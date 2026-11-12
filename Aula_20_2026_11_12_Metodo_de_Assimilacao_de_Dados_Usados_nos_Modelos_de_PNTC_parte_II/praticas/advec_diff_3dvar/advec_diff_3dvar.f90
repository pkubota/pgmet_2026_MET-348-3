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


program advec_diff_3dvar
  use  Class_Fields, only: xa,ua
  USE ModAdvection, only : AnaliticFunction
  USE Class_WritetoGrads, Only : SchemeWriteCtl,SchemeWriteData
  use Model, only: Init,Run_Model,Finalize
  use LinearSolve, only :inverse
  use iso_fortran_env, only: real64
  implicit none
  LOGICAL :: monotone=.true.
  integer, parameter :: nx = 50, nt = 250, nobs=nx
  real(8), parameter :: L = 1.0_real64, dt = 0.001_real64, u = 1.0_real64, K = 0.01_real64
  real(real64) :: obs_error = 0.05_real64
  real(real64) :: model_error = 0.1_real64       ! sigma_b
  real(real64) :: Lcorr = 0.1_real64             ! comprimento de correlação
  INTEGER, PARAMETER  :: fixedSeedValue = 123456 ! A fixed, arbitrary integer value
  real(8) :: step_t,rand1,rand2
  real(8),allocatable :: T(:), Tnew(:), Tb(:), Ta(:), y_obs(:), diff(:),increment(:)
  real(8),allocatable :: B(:,:), Ht(:,:), H(:,:)
  real(8),allocatable :: R1(:,:),R(:,:), Id(:,:), tmp(:,:)
  real(8),allocatable :: mat_A(:,:),Mat_B(:,:),Mat_E(:), gain(:,:),Yobs_HXb(:)
  real(8),allocatable :: x(:)  
  integer, allocatable :: obs_idx(:)
  integer :: i, j,n,idim,jdim,irec,test,n_step,seedSize
  INTEGER, ALLOCATABLE :: seed(:)
  REAL (KIND=8) :: qmax
  REAL (KIND=8) :: qmin

  allocate(x(nx));x=0.0
  allocate(T(nx));T=0.0
  allocate(Tnew(nx));Tnew=0.0
  allocate(Tb(nx));Tb=0.0
  allocate(Ta(nx));Ta=0.0
  allocate(y_obs(nobs));y_obs=0.0
  allocate(diff(nx));diff=0.0
  allocate(B(nx,nx));B=0.0
  allocate(Ht(nx,nobs));Ht=0.0
  allocate(H(nobs,nx));H=0.0
  allocate(R1(nx,nx));R1=0.0
  allocate(R(nobs,nobs));R=0.0
  allocate(Id(nobs,nobs));Id=0.0
  allocate(tmp(nobs,nobs));tmp=0.0
  allocate(mat_A(nx,nx));mat_A=0.0
  allocate(Mat_B(nx,nx));Mat_B=0.0
  allocate(Mat_E(nobs));Mat_E=0.0
  allocate(Yobs_HXb(nobs));Yobs_HXb=0.0
  allocate(gain(nx,nobs));gain=0.0
  allocate(increment(nx));gain=0.0
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
  ! ---------------------------
  ! Construir B (covariancia gaussiana)
  ! B_ij = sigma_b^2 * exp( - (x_i - x_j)^2 / (2 L^2) )
  ! ---------------------------
  B = model_error**1  
  do i = 1, nx
     do j = 1, nx
         call random_number(rand1)
         call random_number(rand2)
         ! Box-Muller ? z ~ N(0,1)
         !
         !     --           --
         !    | b11  b12  b13 |
         !    | b21  b22  b23 | = ! Box-Muller
         !    | b31  b32  b33 |
         !    --            --
         !
         B(i,j) = model_error**1.2 * sqrt((rand1**2) + (rand2**2)) *  exp( - (x(i)-x(j))**2 / (2.0_real64 * Lcorr**2) )
     end do
  end do
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
      R1(obs_idx(i),obs_idx(i)) = 0.02_real64
  END DO
  !
  ! --- Operador de observacao H ---
  !
  H (1:nobs,1:nx)        = 0.0_real64
  do j = 1, size(pack(obs_idx, obs_idx /= -999.0),dim=1)
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
         H(j,obs_idx(i)) = obs_error**1.2 * sqrt((rand1**2) + (rand2**2))*  exp( - (x(i)-x(j))**2 / (2.0_real64 * Lcorr**2) )
     END DO
  END DO
  Ht(1:nx,1:nobs)       = transpose(H(1:nobs,1:nx))
  !
  ! --- Condicao inicial (gaussiana) ---
  !
  irec=0
  step_t=0.0
  n_step=0 
  Tb (1:nx)  = AnaliticFunction(n_step)!c0(1:nx)                 !pluma de concentracao
  T  (1:nx)  = Run_Model(step_t,n_step,real(nt,kind=real64),monotone)!0.8_real64 * c0(1:nx)    !pluma de concentracao
 
  test=SchemeWriteData(T  (1:nx),irec)
  test=SchemeWriteData(Tb (1:nx),irec)
  test=SchemeWriteData(T  (1:nx),irec)

  ! ---------------------------
  ! Simulacao (verdade e modelo)
  ! ---------------------------

  test=0
  !
  ! --- Simulação do modelo (background) ---
  !
  do n = 1, nt
      T (1:nx)  = AnaliticFunction(n_step)                              ! true solution
      Tb(1:nx)  = Run_Model(step_t,n_step,real(nt,kind=real64),monotone)! estado background pluma de concentracao
      ! ---------------------------
      ! Gerar observacoes (ruido Gaussiano) usando Box-Muller
      ! ---------------------------
      call random_seed()
      do i = 1, size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         call random_number(rand1)
         call random_number(rand2)
         if (rand1 <= 0.0_real64) rand1 = MAX(1.0e-12_real64,rand1)
         y_obs(i) = T(obs_idx(i)) + (Tb(obs_idx(i))-T (obs_idx(i)))/2.0
      end do
      !
      ! --- Operador de observacao adaptativo H ---
      !
      H (1:nobs,1:nx)        = 0.0_real64
      DO i=1,size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         do j = 1,  size(pack(obs_idx, obs_idx /= -999.0),dim=1)
            call random_number(rand1)
            call random_number(rand2)
            if( (Tb(obs_idx(i))-T (obs_idx(i)))/2.0 > -0.2*T(obs_idx(i)) .and. &
                (Tb(obs_idx(i))-T (obs_idx(i)))/2.0 <  0.2*T(obs_idx(i))) then  
               H(j,obs_idx(i)) = obs_error * sqrt(((Tb(obs_idx(i))-T (obs_idx(i)))**2)/2)
            endif
         END DO
      END DO
      Ht(1:nx,1:nobs)       = transpose(H(1:nobs,1:nx))
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
      !R (1:nobs,1:nobs)
      !H (1:nobs,1:nobs)
      !Ht(1:nobs,1:nobs)
      !
      !                                                       | 1 2 3 |   | a b |    | x1  X1|
      !                                                       | 4 5 6 | * | c d |  = | x2  x2|
      !                                                       | 7 8 9 |   ! e f |    | x3  x3|
      !                                                              A(nx,nobs)
      !                                                         ---------------------
      !                                                        |     --------        \
      ! --- Aplicação do 3DVar ---                             |    /        \        \
      tmp(1:nobs,1:nobs) = matmul(H(1:nobs,1:nx), matmul(B(1:nx,1:nx),  Ht(1:nx, 1:nobs))) + R(1:nobs,1:nobs)
      !                            \                       A(nx,nobs)
      !                             \                      /
      !                               ---B(nobs,nobs) -----
  
      ! --- Inverse Operator tmp ---
  
      jdim=size(pack(obs_idx, obs_idx /= -999.0),dim=1)
      idim=size(pack(obs_idx, obs_idx /= -999.0),dim=1)
      DO j=1,jdim
         DO i=1,idim
            mat_A(i,j) =tmp (i,j)
            Mat_B(i,j) =tmp (i,j)
         END DO
      END DO

      call inverse(mat_A(1:idim,1:jdim),Mat_B(1:idim,1:jdim),idim)

      DO j=1,jdim
         DO i=1,idim
            Id (i,j)= mat_B(i,j) 
         END DO
      END DO
      ! ---------------------------
      ! Calcular ganho K = B_sub * Sinv   (nx x nobs) * (nobs x nobs) -> (nx x nobs)
      ! ---------------------------
      !                         ------D(nx,nobs)-----------------
      !                                                  ------C(nx,nobs)-----
      gain(1:nx,1:nobs) = matmul(B(1:nx,1:nx), matmul(Ht(1:nx,1:nobs), Id(1:nobs,1:nobs)))   ! ganho de Kalman
      !
      ! atualizar analise
      !
      !                                                       | 1 2 3 |   | a   |    | x1 |
      !                                                       | 4 5 6 | * | c   |  = | x2 |
      !                                                       | 7 8 9 |   ! e   |    | x3 |
      !                                                              A(nx,nobs)
      !
      !
      !                            -----E(nobs,nx)---
      !                           /                  \
      Mat_E(1:nobs)=matmul(H(1:nobs,1:nx) , Tb(1:nx) )
      !
      Yobs_HXb(1:nobs) = (y_obs(1:nobs) - Mat_E(1:nobs))
      !
      !                               --------- F(nx,nx)------
      !                               /                       \
      increment(1:nx) = matmul(gain(1:nx,1:nobs), Yobs_HXb(1:nobs))
      !
      if( monotone )THEN
         Ta(1:nx)  = Tb(1:nx)  + increment(1:nx)
         Tnew(1)=Ta(1)
         qmax = max(Ta(nx),Ta(1))
         qmin = min(Ta(nx),Ta(1))
         Tnew(1) = max(qmin,min(qmax,Ta(1)))
         do i = 2, nx
            qmax = max(Ta(i-1),Ta(i))
            qmin = min(Ta(i-1),Ta(i))
            Tnew(i) = max(qmin,min(qmax,Ta(i)))
         end do
         Ta=Tnew
      else
        Ta(1:nx)  = Tb(1:nx)  + increment(1:nx)
      end if

      test=SchemeWriteData(Tb          (1:nx),irec)
      test=SchemeWriteData(T           (1:nx),irec)
      test=SchemeWriteData(Ta          (1:nx),irec)
      !
      ! --- Resultado ---
      !
      DO i=1,size(pack(obs_idx, obs_idx /= -999.0),dim=1)
         print *, "----------------"
         print *, "Ponto observado:", obs_idx(i)
         print *, "Background:     ", Tb(obs_idx(i))
         print *, "Observacao:     ", y_obs(i)
         print *, "Analise:        ", Ta(obs_idx(i))
      END DO
  end do
  test=SchemeWriteCtl(n_step)


contains

  subroutine advec_diff_step(T, Tnew, u, K, dx, dt, nx)
    implicit none
    integer, intent(in) :: nx
    real(8), intent(in) :: u, K, dx, dt
    real(8), dimension(nx), intent(in) :: T
    real(8), dimension(nx), intent(out) :: Tnew
    integer :: i
    Tnew = T
    do i = 2, nx-1
       Tnew(i) = T(i) - (u*dt/(2._real64*dx))*(T(i+1)-T(i-1)) + (K*dt/(dx*dx))*(T(i+1)-2*T(i)+T(i-1))
    end do
    ! Condicaes de contorno periodicas
    Tnew(1) = T(nx-1)
    Tnew(nx) = T(2)
  end subroutine advec_diff_step

  function inv1x1(a) result(ainv)
    real(8), dimension(1,1), intent(in) :: a
    real(8), dimension(1,1) :: ainv
    ainv(1,1) = 1.0_real64 / a(1,1)
  end function inv1x1

end program advec_diff_3dvar
