MODULE LinearSolve 
 IMPLICIT NONE
 PRIVATE
 INTEGER, PARAMETER :: r8 = selected_real_kind(15, 307)
 INTEGER, PARAMETER :: r4 = selected_real_kind(6, 37)
 PUBLIC :: solve_tridiag,inverse
CONTAINS
subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
!       a - sub-diagonal (diagonal abaixo da diagonal principal)
!       b - diagonal principal
!       c - sup-diagonal (diagonal acima da diagonal principal)
!       d - parte à direita
!       x - resposta
!       n - número de equações
        integer, intent(in) :: n
        real (r8), dimension (n),intent (in   ) :: a,b,c,d
        real (r8), dimension (n),intent (out) :: x
        real (r8), dimension (n) :: cp, dp
        real (r8) :: m
        integer :: i
! inicializar c-primo e d-primo
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! resolver para vetores c-primo e d-primo
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! inicializar x
         x(n) = dp(n)
! resolver para x a partir de vetores c-primo e d-primo
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
    end subroutine solve_tridiag

  subroutine inverse(a_in,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real (r8) a_in(n,n), c(n,n)
real (r8) L(n,n), U(n,n), b(n), d(n), x(n)
real (r8) coeff
integer i, j, k
real (r8) a(n,n)

a=a_in
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse



END MODULE LinearSolve



! kf_advec_diff.f90
program kf_advec_diff
  use  Class_Fields        , only: xa,ua
  !use  Class_Fields_Adjoint, only: xa,ua
  USE ModAdvection        , only : AnaliticFunction
  USE ModAdvection_Adjoint, only : AnaliticFunction_Adjoint
  USE Class_WritetoGrads, Only : SchemeWriteCtl,SchemeWriteData
  use Model        , only: Init,Run_Model,Finalize
  use Model_Adjoint, only: Init_Adjoint,Run_Model_Adjoint,Finalize_Adjoint
  use LinearSolve, only :inverse
  use iso_fortran_env, only: real64
  implicit none
  integer , parameter :: nx = 60, nt = 400, nobs=nx
  real(real64), parameter :: L = 1.0_real64, u = 1.0_real64, kappa = 0.01_real64
  real(real64), parameter :: dx = L/(nx-1), dt = 0.25_real64*dx
  real(real64), parameter :: sigma_b = 0.10_real64
  real(real64), parameter :: corr_len = 0.15_real64
  real(real64), parameter :: sigma_obs = 0.02_real64
  real(real64), parameter :: q_diag = 1.0e-5_real64
  real(real64) :: Lcorr = 0.1_real64             ! comprimento de correlação
  LOGICAL :: monotone=.true.
  !real(real64)            :: obs_error = 0.01_real64
  !real(real64), parameter :: model_error = 0.1_real64       ! sigma_b
  real(8),allocatable :: x(:)
  integer, allocatable :: obs_idx(:)
  real(real64), allocatable :: A(:,:), truth(:), forecast(:), analysis(:)
  real(real64), allocatable :: x_new(:), x_true(:), x_f(:),mat_b(:,:),mat_a(:,:),mat_invb(:,:), x_a(:), x_b(:)
  real(real64), allocatable :: P(:,:), Pf(:,:), Qa(:,:), Q(:,:),R1(:,:), R(:,:), Rinv(:,:),mat_aux(:,:)
  real(real64), allocatable :: H(:,:), Ht(:,:)
  real(real64), allocatable :: B(:,:), Binv(:,:), mat_Id(:,:)
  real(real64), allocatable :: y_obs(:)
  real(real64), allocatable :: S(:,:), Sinv(:,:)
  real(real64), allocatable :: K(:,:), innov(:)
  integer :: i,j,kk,t, outunit, rmseunit,seedSize,n_step,irec,test
  INTEGER, PARAMETER  :: fixedSeedValue = 123456 ! A fixed, arbitrary integer value
  INTEGER, ALLOCATABLE :: seed(:)
  real(real64) :: step_t,rand1,rand2,d,qmax,qmin

  !call seed_rng(11111)
  allocate(x(nx));x=0.0
  allocate(A(nx,nx))
  allocate(truth(nx),x_new(nx), forecast(nx), analysis(nx))
  allocate(x_true(nx), x_f(nx), x_a(nx),mat_b(nx,nx),mat_a(nx,nx),mat_invb(nx,nx),mat_aux(nx,nx),x_b(nx))
  allocate(P(nx,nx), Pf(nx,nx), Q(nx,nx),R1(nx,nx), R(nobs,nobs), Rinv(nobs,nobs))
  allocate(H(nobs,nx), Ht(nx,nobs))
  allocate(S(size(R,dim=1),size(R,dim=2)), Sinv(size(R,dim=1),size(R,dim=2))   )
  allocate(B(nx,nx), Binv(nx,nx),mat_Id(nx,nx))
  allocate(y_obs(nobs))
  allocate(K(nx,size(R,dim=1)), innov(size(R,dim=1)))
  allocate(obs_idx(nx))
 
  !allocate(obs_idx(nobs))
  ! Inquire about the size of the seed array required by the system
  CALL RANDOM_SEED (size=seedSize)
  ! Allocate the seed array
  ALLOCATE(seed(seedSize))
  ! Initialize the seed array with a consistent value.
  ! This example fills the array with the same fixed value.
  ! More complex initialization based on the fixedSeedValue could also be used.
  seed = fixedSeedValue
  ! Set the random number generator's seed
  CALL RANDOM_SEED(PUT=seed)

  CALL Init(nx)
  CALL Init_Adjoint(nx)
  ! ---------------------------
  ! Grade e condicao inicial
  ! ---------------------------
  do i = 1, nx
     x(i) = xa(i)   ! USE Class_Fields, only : xa real(i-1, real64) / real(nx-1, real64)
  end do
  !
  ! Ponto observado: value of the observation
  !
  obs_idx=-999
  j=1
  do i = 1, nx,1
      obs_idx(i) = MIN(j,nx)   ! Fortran 1-based
      j = j + REAL(nx)/REAL(nobs)
      if(j>nx)exit
  end do

  if(size(pack(obs_idx, obs_idx /= -999.0),dim=1) /= nobs) then
      stop 'error'
  end if 

  call build_A(A, nx, dx, dt, u, kappa) 

   !
   !call build_B_exponential(B, nx, dx, sigma_b, corr_len)
   ! ---------------------------
   ! Construir B (covariancia gaussiana)
   ! B_ij = sigma_b^2 * exp( - (x_i - x_j)^2 / (2 L^2) )
   ! ---------------------------
   B = sigma_b**1  
   do i = 1, nx
      do j = 1, nx
          call random_number(rand1)
          call random_number(rand2)
          ! Box-Muller ? z ~ N(0,1)
          !        B(i,j) = sigma_b**2 * exp(-d/corr_len)
          !     --           --
          !    | b11  b12  b13 |
          !    | b21  b22  b23 | = ! Box-Muller
          !    | b31  b32  b33 |
          !    --            --
          !
          !B(i,j) = model_error**1.2 * sqrt((rand1**2) + (rand2**2)) *  exp( - (x(i)-x(j))**2 / (2.0_real64 * Lcorr**2) )
          d = x(i) - x(j)
          B(i,j) = sigma_b**2 * exp( - (d*d) / (2.0_real64 * Lcorr**2) )
      end do
   end do  
   call invert_matrix(B, Binv, nx)

   ! initial P = B
   !
   !     --           --
   !    | b11  b12  b13 |
   !    | b21  b22  b23 | = ! Box-Muller
   !    | b31  b32  b33 |
   !    --            --
   !
   ! Construir B (covariancia gaussiana)
   P = B
   !
   ! Q (model error) small diagonal
   !
   
   Q = 0.0_real64
   do i = 1, nx
     Q(i,i) = q_diag   !q_diag = 1.0e-5_real64
   end do
  !
  ! --- Operador de observacao H ---
  ! H(i,j) = 1 se o i-ésimo ponto de observação é o j-ésimo ponto da grade
  ! H: nobs x nx
  !
  H (1:nobs,1:nx)        = 0.0_real64
  !do j = 1, size(pack(obs_idx, obs_idx /= -999.0),dim=1)
     DO i = 1, size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         call random_number(rand1)
         call random_number(rand2)
         !
         !     --           --
         !    | R11  R12  R13 |
         !    | R21  R22  R23 | = ! Box-Muller
         !    | R31  R32  R33 |
         !    --            --
         !
         !
         !              grid model
         !            --           --
         !           | R11  R12  R13 |
         ! pts obs   | R21  R22  R23 | = ! Box-Muller
         !           | R31  R32  R33 |
         !           --            --
         !
         H(i, obs_idx(i)) = 1.0_real64
     END DO
  !END DO
  Ht(1:nx,1:nobs)       = transpose(H(1:nobs,1:nx))
  !
  ! --- Matriz de erro de observacao R ---  
  !
  R1(1:nx,1:nx) = 0.0_real64
  DO i=1,size(pack(obs_idx, obs_idx /= -999.0),dim=1)
      !
      !              grid model
      !            --           --
      !           | R11   0    0  |
      ! pts obs   |  0   R22   0  | = ! Box-Muller
      !           |  0    0   R33 |
      !           --            --
      !
      !R1(i,i) = sigma_obs**2

      R1(obs_idx(i),obs_idx(i)) = sigma_obs**2
  END DO
  !
  ! --- Condicao inicial (gaussiana) ---
  !
   irec=0
   step_t=0.0
   n_step=0 
   truth (1:nx)  = AnaliticFunction(n_step)!c0(1:nx)                 !pluma de concentracao
   x_a(1:nx)     =truth (1:nx) !Run_Model(step_t,n_step,real(nt,kind=real64),monotone)!0.8_real64 * c0(1:nx)    !pluma de concentracao
   x_b         = x_a
   forecast(:) = x_a
   analysis(:) = x_a

   !truth(1,:) = x_true
   !do t = 1, nt
   !  truth(:) = matmul(A, truth(:))
   !end do
   test=SchemeWriteData(forecast  (1:nx),irec)
   test=SchemeWriteData(truth     (1:nx),irec)
   test=SchemeWriteData(analysis  (1:nx),irec)
   
   ! ---------------------------
   ! Simulacao (verdade e modelo)
   ! ---------------------------

   test=0
   !
   ! --- Simulação do modelo (background) ---
   !

   ! Loop time: forecast + update
   do t = 1, nt
      !
      ! forecast
      !
      x_f           = Run_Model(step_t,n_step,real(nt,kind=real64),monotone)!0.8_real64 * c0(1:nx)    !pluma de concentracao
      truth (1:nx)  = AnaliticFunction(n_step)!c0(1:nx)                 !pluma de concentracao
      !
      ! Pf  Matrix da covarinacia do erro da previsao
      !
      Pf = matmul(A(1:nx,1:nx), matmul(P(1:nx,1:nx), transpose(A(1:nx,1:nx)))) + Q(1:nx,1:nx)

      forecast(:) = x_f
      ! ---------------------------
      ! Gerar observacoes (ruido Gaussiano)   ! observations usando Box-Muller
      ! ---------------------------
      !call random_seed()
      do i = 1, size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         call random_number(rand1)
         call random_number(rand2)
         if (rand1 <= 0.0_real64) rand1 = MAX(1.0e-12_real64,rand1)
         y_obs(i) = truth(obs_idx(i)) + sigma_obs *((x_f(obs_idx(i))-truth (obs_idx(i)))/2.0)
      end do
      !
      !
      ! --- copia da Matriz de erro de observacao R1 ---  
      !
      !
      R(1:nobs,1:nobs) = 0.0_real64
      DO j=1,size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         DO i=1,size(pack(obs_idx, obs_idx /= -999.0),dim=1)
          !              grid model
          !            --           --
          !           | R11   0    0  |
          ! pts obs   |  0   R22   0  | = ! Box-Muller
          !           |  0    0   R33 |
          !           --            --
          !
          R(i,j)= R1(obs_idx(i),obs_idx(j)) 
        END DO
      END DO
      Rinv(1:nobs,1:nobs)       = transpose(R(1:nobs,1:nobs))

      ! Kalman gain
      ! S = H Pf H^T + R            =>  (m x m)
      !
      S(1:nobs,1:nobs) = matmul(H(1:nobs,1:nx), matmul(Pf(1:nx,1:nx), transpose(H(1:nobs,1:nx)))) + R(1:nobs,1:nobs)
      call invert_matrix(S(1:nobs,1:nobs), Sinv(1:nobs,1:nobs), size(R,1))
      !
      ! K = Pf H^T S^{-1}            => (nx x m)
      !
      K(1:nx,1:nobs) = matmul(Pf(1:nx,1:nx), matmul(transpose(H(1:nobs,1:nx)), Sinv(1:nobs,1:nobs)))
      !
      innov(1:nobs) = y_obs(1:nobs) - matmul(H(1:nobs,1:nx), x_f(1:nx))

      
      if( monotone )THEN
         x_a(1:nx) = x_f(1:nx) + matmul(K(1:nx,1:nobs), innov(1:nobs))
         x_new(1)=x_a(1)
         qmax = max(x_a(nx),x_a(1))
         qmin = min(x_a(nx),x_a(1))
         x_new(1) = max(qmin,min(qmax,x_a(1)))
         do i = 2, nx
            qmax = max(x_a(i-1),x_a(i))
            qmin = min(x_a(i-1),x_a(i))
            x_new(i) = max(qmin,min(qmax,x_a(i)))
         end do
         x_a=x_new
      else
        x_a(1:nx) = x_f(1:nx) + matmul(K(1:nx,1:nobs), innov(1:nobs))
      end if
      
      ! I = identidade(nx)
      
      mat_aux(1:nx,1:nx) = identity(nx) - matmul(K(1:nx,1:nobs), H(1:nobs,1:nx))

      P(1:nx,1:nx) = matmul( mat_aux(1:nx,1:nx), matmul(Pf(1:nx,1:nx), transpose(mat_aux(1:nx,1:nx))) ) &
                   + matmul( K(1:nx,1:nobs), matmul(R(1:nobs,1:nobs), transpose(K(1:nx,1:nobs))) )

      analysis(1:nx) = x_a(1:nx)
      
      !P(1:nx,1:nx) = matmul( (identity(nx) - matmul(K(1:nx,1:nobs),H(1:nobs,1:nx))), Pf(1:nx,1:nx) )
      !analysis(1:nx) = x_a(1:nx)


      print*,' ---- '
      write(*,'(I6,1X,60F12.6,1X,60F12.6,1X,60F12.6)') t-1, truth(:), forecast(:), analysis(:)
      test=SchemeWriteData(forecast          (1:nx),irec)
      test=SchemeWriteData(truth             (1:nx),irec)
      test=SchemeWriteData(analysis          (1:nx),irec)
 
  end do
  test=SchemeWriteCtl(n_step)


contains
  ! ---------------- truth IC ----------------
  subroutine make_truth_ic(nx, dx, x0)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dx
    real(real64), intent(out) :: x0(nx)
    integer :: i
    do i = 1, nx
      x0(i) = exp( - (( (real(i-1,real64)*dx) - 0.25_real64 )**2) / (2.0_real64*(0.05_real64**2)) ) + &
              0.5_real64 * exp( - (( (real(i-1,real64)*dx) - 0.65_real64 )**2) / (2.0_real64*(0.08_real64**2)) )
    end do
  end subroutine make_truth_ic

  subroutine build_A(A, nx, dx, dt, u, kappa)

    integer, intent(in) :: nx
    real(real64), intent(in) :: dx, dt, u, kappa
    real(real64), intent(out) :: A(nx,nx)
    integer :: i
    !A(1,1)   = 1.0_real64
    !A(nx,nx) = 1.0_real64
    !
    ! A(n+1,i) - A(n,i) = - (udt/dx)*A(n,i+1) - (udt/dx)*A(n,i)  + (k*dt/dxdx)*A(n,i+1) + (k*dt/dxdx)*2A(n,i) + (k*dt/dxdx)*A(n,i-1)
    !-----------------
    !
    !
    ! A(n+1,i) = A(n,i) - (udt/dx)*A(n,i+1) - (udt/dx)*A(n,i)  + (k*dt/dxdx)*A(n,i+1) + (k*dt/dxdx)*2A(n,i) + (k*dt/dxdx)*A(n,i-1)
    !  
    ! A(n+1,i) = A(n,i) - (udt/dx)*A(n,i+1) - (udt/dx)*A(n,i)  + (k*dt/dxdx)*A(n,i+1) + (k*dt/dxdx)*2A(n,i) + (k*dt/dxdx)*A(n,i-1)
    !
    ! A(n+1,i) = (k*dt/dxdx)*A(n,i-1) + A(n,i) - (udt/dx)*2A(n,i) + (k*dt/dxdx)*A(n,i) - (udt/dx)*A(n,i+1)  + (k*dt/dxdx)*A(n,i+1)
    !
    ! A(n+1,i) = (k*dt/dxdx)*A(n,i-1) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,i) - [(udt/dx) + (k*dt/dxdx) ] * A(n,i+1)
    !
    ! k*dt/dxdx)*A(n,i-1) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,i) - [(udt/dx) + (k*dt/dxdx) ] * A(n,i+1) = A(n+1,i)
    !
    !(i= 1) =>  k*dt/dxdx)*A(n,0) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,1) - [(udt/dx) + (k*dt/dxdx) ] * A(n,2) = A(n+1,1)
    !(i= 2) =>  k*dt/dxdx)*A(n,1) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,2) - [(udt/dx) + (k*dt/dxdx) ] * A(n,3) = A(n+1,2)
    !(i= 3) =>  k*dt/dxdx)*A(n,2) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,3) - [(udt/dx) + (k*dt/dxdx) ] * A(n,4) = A(n+1,3)
    !(i= 4) =>  k*dt/dxdx)*A(n,3) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,4) - [(udt/dx) + (k*dt/dxdx) ] * A(n,5) = A(n+1,4)
    !
    !(i= i) =>  k*dt/dxdx)*A(n,i-1) + [1 - (udt/dx) + (2k*dt/dxdx)]*A(n,i) - [(udt/dx) + (k*dt/dxdx) ] * A(n,i+1) = A(n+1,i)
    !
    ! fronteria nao ciclica
    !   --                                                                             --         --    --      --      --
    !  | k*dt/dxdx)      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)  0  0   |    |A(n,0  ) |   | A(n+1,1) |
    !  | k*dt/dxdx)      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)         |    |A(n,1  ) |   | A(n+1,2) |
    !  | k*dt/dxdx)      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)         |  * |A(n,2  ) | = | A(n+1,3) |
    !  | k*dt/dxdx)      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)         |    |A(n,3  ) |   | A(n+1,4) |
    !
    !  | k*dt/dxdx)      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)         |    |A(n,i-1) |   | A(n+1,i) |
    !   --                                                                             --         --     --     --      --
    !
    ! fronteria ciclica
    !   --                                                                                                                                             --     --    --      --      --
    !  | [1 - (udt/dx) + (2k*dt/dxdx)]            [(udt/dx) + (k*dt/dxdx)]                    0                       [k*dt/dxdx)]                      |    |A(n,0  ) |   | A(n+1,1) |
    !  | [ k*dt/dxdx)]                            [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)]          0                                 |    |A(n,1  ) |   | A(n+1,2) |
    !  |      0                                   [k*dt/dxdx)]                      [1 - (udt/dx) + (2k*dt/dxdx)]     [(udt/dx) + (k*dt/dxdx)]          |  * |A(n,2  ) | = | A(n+1,3) |
    !  |      0                                   0                                 [k*dt/dxdx) ]                     [1 - (udt/dx) + (2k*dt/dxdx)]     |    |A(n,3  ) |   | A(n+1,4) |
    !
    !  | [(udt/dx) + (k*dt/dxdx)]                 0                                 k*dt/dxdx)                        [1 - (udt/dx) + (2k*dt/dxdx)]     |    |A(n,i-1) |   | A(n+1,i) |
    !   --                                                                                                                                            
    A(1,2  ) = u*dt/dx + kappa*dt/(dx*dx)
    A(1,1  ) = 1.0_real64 - u*dt/dx - 2.0_real64*kappa*dt/(dx*dx)
    !A(1,0) = 1.0_real64
    do i = 2, nx-1
      A(i,i+1) = u*dt/dx + kappa*dt/(dx*dx)
      A(i,i)   = 1.0_real64 - u*dt/dx - 2.0_real64*kappa*dt/(dx*dx)
      A(i,i-1) = kappa*dt/(dx*dx)
    end do


    !A(nx,nx+1  ) = 0
    A(nx,nx    ) = 1.0_real64 - u*dt/dx - 2.0_real64*kappa*dt/(dx*dx)
    A(nx,nx-1  ) =  u*dt/dx + kappa*dt/(dx*dx)

  end subroutine build_A

  subroutine build_B_exponential(B, nx, dx, sigma_b, corr_len)
    real(real64), intent(out) :: B(nx,nx)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dx, sigma_b, corr_len
    integer :: i,j
    real(real64) :: xi, xj, d
    do i = 1, nx
      xi = (real(i-1,real64) * dx)
      do j = 1, nx
        xj = (real(j-1,real64) * dx)
        d = abs(xi-xj)
        B(i,j) = sigma_b**2 * exp(-d/corr_len)
      end do
    end do
  end subroutine build_B_exponential

  subroutine compute_K_and_update(Pf, H, Ht, R, Rinv, xf, y_obs, xa, Pa)
    real(real64), intent(in) :: Pf(nx,nx), H(:,:), Ht(:,:), R(:,:), Rinv(:,:)
    real(real64), intent(in) :: xf(nx), y_obs(:)
    real(real64), intent(inout) :: xa(nx), Pa(nx,nx)
    real(real64) :: S(size(R,1),size(R,2)), Sinv(size(R,1),size(R,2))
    real(real64) :: K(nx,size(R,1)), innov(size(R,1))
    integer :: info
    S = matmul(H, matmul(Pf, transpose(H))) + R
    call invert_matrix(S, Sinv, size(R,1))
    K = matmul(Pf, matmul(transpose(H), Sinv))
    innov = y_obs - matmul(H, xf)
    xa = xf + matmul(K, innov)
    Pa = matmul( (identity(nx) - matmul(K,H)), Pf )
  end subroutine compute_K_and_update

  function identity(n) result(Id)
    integer, intent(in) :: n
    real(real64) :: Id(n,n)
    integer :: ii
    Id = 0.0_real64
    do ii=1,n; Id(ii,ii)=1.0_real64; end do
  end function identity
  ! ---------------- matrix inversion (Gauss-Jordan pivoting) ----------------

  subroutine invert_matrix(A, Ainv, n)
    real(real64), intent(in) :: A(n,n)
    real(real64), intent(out) :: Ainv(n,n)
    integer, intent(in) :: n
    real(real64), allocatable :: aug(:,:)
    integer :: i,j,k,piv
    real(real64) :: maxv,fac,tmp
    allocate(aug(n,2*n))
    aug(:,1:n) = A
    aug(:,n+1:2*n) = 0.0_real64
    do i=1,n; aug(i,n+i)=1.0_real64; end do
    do i=1,n
      ! pivot selection (partial pivoting)
      piv = i; maxv = abs(aug(i,i))
      do k=i+1,n
        if (abs(aug(k,i)) > maxv) then; maxv = abs(aug(k,i)); piv = k; end if
      end do
      if (piv /= i) then
        aug([i,piv],:) = aug([piv,i],:)
      end if
      if (abs(aug(i,i)) < 1.0e-16_real64) aug(i,i) = 1.0e-16_real64
      ! normalize row i
      fac = aug(i,i)
      aug(i,:) = aug(i,:)/fac
      ! eliminate others
      do j=1,n
        if (j /= i) then
          tmp = aug(j,i); aug(j,:) = aug(j,:) - tmp*aug(i,:)
        end if
      end do
    end do
    Ainv = aug(:, n+1:2*n)
    deallocate(aug)
  end subroutine invert_matrix

  function randn() result(r)
    real(real64) :: r, u1,u2
    call random_number(u1); call random_number(u2)
    if (u1 < 1.0e-16_real64) u1 = 1.0e-16_real64
    r = sqrt(-2.0_real64*log(u1))*cos(2.0_real64*acos(-1.0_real64)*u2)
  end function randn

  function randn_vec(n) result(v)
    integer, intent(in):: n
    real(real64) :: v(n); integer ii
    do ii=1,n; v(ii) = randn(); end do
  end function randn_vec

  subroutine seed_rng(seed)
    integer, intent(in):: seed
    integer :: n,i; integer, allocatable :: s(:)
    call random_seed(size=n); allocate(s(n))
    do i=1,n; s(i)=i+seed; end do
    call random_seed(put=s); deallocate(s)
  end subroutine seed_rng

  function rmse_vec(a,b) result(rmse)
    real(real64), intent(in):: a(:), b(:); real(real64) :: rmse; integer n
    n = size(a); rmse = sqrt(sum((a-b)**2)/real(n,real64))
  end function rmse_vec

end program kf_advec_diff
