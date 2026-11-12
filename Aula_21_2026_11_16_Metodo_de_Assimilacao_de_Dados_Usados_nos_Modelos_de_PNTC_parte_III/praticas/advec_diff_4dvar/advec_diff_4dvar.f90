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
! advec_diff_4dvar.f90
! 4D-Var (strong-constraint) for 1D advection-diffusion with:
!  - Exponential B (non-diagonal)
!  - Preconditioned Conjugate Gradient (Jacobi preconditioner using diag(LHS))
!  - Configurable window length (win) and stride (stride)
!
program advec_diff_4dvar
  use LinearSolve, only :inverse
  implicit none
  integer , parameter :: dp = selected_real_kind(15,307)
  integer , parameter :: nx = 60            ! grid points
  integer , parameter :: nt = 120           ! total time steps
  integer , parameter :: win = 12           ! window length (steps)
  real(dp), parameter :: L = 1.0_dp
  real(dp), parameter :: u = 1.0_dp
  real(dp), parameter :: kappa = 0.01_dp
  real(dp), parameter :: dx = L / real(nx-1, dp)
  real(dp), parameter :: dt = 0.25_dp * dx   ! choose stable dt
  real(dp), parameter :: sigma_b = 0.10_dp
  real(dp), parameter :: sigma_obs = 0.02_dp
  INTEGER, PARAMETER  :: fixedSeedValue = 123456 ! A fixed, arbitrary integer value

  ! Observations: indices where observations exist (m per time)
  integer, parameter :: m = 6
  integer, parameter :: obs_idx(m) = (/ 6, 14, 22, 30, 38, 46 /)!nx = 60            ! grid points

  ! arrays
  real(dp), allocatable :: A(:,:)    , Ap(:,:,:), AT(:,:)      , ATp(:,:,:)
  real(dp), allocatable :: truth(:,:), bkg(:,:) , analysis(:,:)
  real(dp), allocatable :: x0_true(:), x0_b(:), x0_rand(:)
  real(dp), allocatable :: y(:,:)
  real(dp), allocatable :: B(:,:), Binv(:,:), R(:,:), Rinv(:,:),mat_A(:,:),Mat_B(:,:)
  real(dp), allocatable :: H(:,:), Ht(:,:)
  real    , allocatable :: dump_var(:)
  INTEGER, ALLOCATABLE :: seed(:)

  ! temp / CG matrices
  real(dp), allocatable :: LHS(:,:), RHS(:), x0_opt(:)

  integer :: i,j,k,t,win_k,info,jdim,idim,seedSize,lrec,irec,icount
  integer :: outunit,outunitbin, rmseunit
  real(dp) :: t0,num
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

  ! allocate
  allocate(A (nx,nx)      , AT (nx,nx))
  allocate(Ap(nx,nx,win+1), ATp(nx,nx,win+1))
  allocate(truth(nt+1,nx) , bkg(nt+1,nx), analysis(nt+1,nx))
  allocate(x0_true(nx)    , x0_b(nx),x0_rand(nx))
  allocate(y(nt,m))
  allocate(B(nx,nx)       , Binv(nx,nx) , R(m,m), Rinv(m,m))
  allocate(mat_A(nx,nx)   , Mat_B(nx,nx))
  allocate(H(m,nx)        , Ht(nx,m))
  allocate(LHS(nx,nx)     , RHS(nx)     , x0_opt(nx),dump_var(nx))
  !
  ! build A
  !
  call build_A(A, nx, dx, dt, u, kappa)
  !
  !  --           --
  ! | d   up   0   0|
  ! | dw  d   up   0|   =  I = > Ap
  ! | 0   dw   d  up|
  ! | 0   0   dw   d|
  !  --           --  
  !
  AT = transpose(A)
  !
  ! build powers A^k and (A^T)^k
  !
  !  --     --
  ! | 1  0  0 |
  ! | 0  1  0 |   =  I = > Ap
  ! | 0  0  1 |
  !  --     --
  Ap (:,:,1)  = identity(nx)
  !
  !  --     --
  ! | 1  0  0 |
  ! | 0  1  0 |   =  I = > ATp
  ! | 0  0  1 |
  !  --     --
  !
  ATp(:,:,1)  = identity(nx)
  !
  
  do k = 2, win + 1 !win window length (steps)
      !  --     --     --          --
      ! | 1  0  0 |   | d   up   0  |      | d   up    0|
      ! | 0  1  0 | * | dw  d   up  |   =  | dw  d    up|  = > Ap
      ! | 0  0  1 |   | 0   dw   d  |      | 0   dw    d|
      !  --     --     --         --

      Ap (:,:,k)  = matmul   (A, Ap(:,:,k-1))

      ATp(:,:,k)  = transpose(Ap(:,:,k))

  end do

  ! truth initial condition (two gaussian bumps)

  call make_truth_ic(nx, dx, x0_true)

  truth(1 ,:) = x0_true(1:nx)

  do t = 1, nt
    truth(t+1,:) = matmul(A, truth(t,:))
  end do

  ! observations (one per time step at obs indices)

  do t = 1, nt
    do i = 1, m ! Observations: indices where observations exist (m per time)
       CALL RANDOM_NUMBER(num)
       y(t,i) = truth(t+1, obs_idx(i)) + sigma_obs *num ! randn()
    end do
  end do

  ! background: perturb initial truth and propagate free

  x0_b(1:nx) = x0_true(1:nx) + sigma_b * randn_vec(nx)

  t=1 
  bkg(t,1:nx) = x0_b(1:nx)

  do t = 1, nt
    bkg(t+1,1:nx) = matmul(A, bkg(t,1:nx))
  end do

  ! build B (diagonal) and Binv

  B = 0.0_dp
  Binv = 0.0_dp
  do i = 1, nx
    B(i,i)    = sigma_b**2
  end do
  ! --- Inverse Operator tmp ---
  
  mat_A= 0.0_dp
  Mat_B= 0.0_dp
  jdim=size(B,dim=2)
  idim=size(B,dim=1)
  DO i=1,idim
      mat_A(i,i) = B(i,i) 
      Mat_B(i,i) = 0.0_dp
  END DO
  call inverse(mat_A(1:idim,1:jdim),Mat_B(1:idim,1:jdim),idim)
  DO j=1,jdim
     DO i=1,idim
         Binv(i,j)= mat_B(i,j) 
     END DO
  END DO
  ! build R and Rinv (m x m)
  R = 0.0_dp
  Rinv = 0.0_dp
  do i = 1, m
      R   (i,i)    = sigma_obs**2 
  end do
  ! --- Inverse Operator tmp ---
  
  mat_A= 0.0_dp
  Mat_B= 0.0_dp
  jdim=size(R,dim=2)
  idim=size(R,dim=1)
  DO i=1,idim
      mat_A(i,i) = R(i,i) 
      Mat_B(i,i) = 0.0_dp
  END DO
  call inverse(mat_A(1:idim,1:jdim),Mat_B(1:idim,1:jdim),idim)
  DO j=1,jdim
     DO i=1,idim
       Rinv(i,j)= mat_B(i,j) 
     END DO
  END DO

  ! build H (m x nx) and Ht (nx x m)
  H = 0.0_dp
  do i = 1, m
    H(i, obs_idx(i)) = 1.0_dp
  end do
  Ht = transpose(H)

  ! open output files 
  INQUIRE(iolength=lrec)dump_var 
  open(newunit=outunitbin, file="4dvar_out.bin",access='direct',recl=lrec,form='unformatted',&
       action="write", status="unknown")
  open(newunit=outunit, file="4dvar_out.dat", status="replace", action="write")
  write(outunit, '(A)') '# t truth(...) background(...) analysis(...)'
  open(newunit=rmseunit, file="rmse_4dvar.dat", status="replace", action="write")
  write(rmseunit, '(A)') '# t rmse_background rmse_analysis'
  open(newunit=outunit, file="4dvar_out.dat", status="replace", action="write")

  ! initialize analysis with background initial guess
  analysis(1,1:nx) = bkg(1,1:nx)
  irec=0
  icount=0
  t = 1
  do while (t <= nt)
    ! Determine window end
    win_k = min(win, nt - t + 1)

    ! Build LHS and RHS for window starting at time t
    ! LHS = Binv + sum_{k=1..win_k} AT^k * H^T * Rinv * H * A^k
    LHS(1:nx,1:nx) = Binv(1:nx,1:nx)
    do k = 1, win_k
      LHS(1:nx,1:nx) = LHS(1:nx,1:nx) + matmul( ATp(1:nx,1:nx,k+1), matmul( matmul(Ht(1:nx,1:m), Rinv(1:m,1:m)), &
                                                                            matmul(H (1:m,1:nx), Ap  (1:nx,1:nx,k+1)) ) )
    end do

    ! RHS = Binv * xb0 + sum_k AT^k * H^T * Rinv * y_k
    RHS(1:nx) = 0.0_dp
    RHS(1:nx) = matmul(Binv(1:nx,1:nx), bkg(t,1:nx))   ! xb0 = background at window start
    do k = 1, win_k
      RHS(1:nx) = RHS(1:nx) + matmul( ATp(1:nx,1:nx,k+1), matmul(Ht(1:nx,1:m), &
                              matmul( Rinv(1:m,1:m)     , y(t + k - 1, 1:m))) )
    end do

    ! Solve LHS * x0_opt = RHS using CG (x0 initial guess = xb0)
    x0_opt(1:nx) = bkg(t,1:nx)    ! initial guess
    
    call cg_solve(LHS(1:nx,1:nx), RHS(1:nx), x0_opt(1:nx), 1000, 1.0e-9_dp)

    ! propagate optimal x0 through the window and store analysis
    do k = 0, win_k-1
      analysis(t + k, 1:nx) = matmul(Ap(1:nx,1:nx,k+1), x0_opt(1:nx))
    end do
    ! also write the last time of window
    analysis(t + win_k, 1:nx) = matmul(Ap(1:nx,1:nx,win_k+1), x0_opt(1:nx))

    ! advance by window length (no overlap here, could do overlap if desired)
    t = t + win_k
    icount = icount + 1
    do i = 1, m
       write(*,  '(I6,1X,60F12.6,1X,60F12.6,1X,60F12.6)') t-1, truth(t,obs_idx(i)), bkg(t,obs_idx(i)), analysis(t,obs_idx(i))
    end do
    print*,'t',icount
    irec=irec+1
    dump_var(1:nx)=  truth(t,1:nx)
    WRITE(outunitbin,rec=irec) dump_var
    irec=irec+1
    dump_var(1:nx)=   bkg(t,1:nx)
    WRITE(outunitbin,rec=irec) dump_var
    irec=irec+1
    dump_var(1:nx)= analysis(t,1:nx)
    WRITE(outunitbin,rec=irec) dump_var

  end do
  Close(outunitbin) 

  ! compute rmse per time and write outputs
  do t = 1, nt+1
    write(outunit,  '(I6,1X,60F12.6,1X,60F12.6,1X,60F12.6)') t-1, truth(t,:), bkg(t,:), analysis(t,:)
    write(rmseunit, '(I6,2F12.6)') t-1, rmse_vec(bkg(t,:), truth(t,:)), rmse_vec(analysis(t,:), truth(t,:))
  end do
  close(outunit)
  close(rmseunit)
  Close(outunitbin) 
  print *, '4D-Var run finished. Output: 4dvar_out.dat, rmse_4dvar.dat'

contains

  !============================================================
  subroutine build_A(A, nx, dx, dt, u, kappa)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt, u, kappa
    real(dp), intent(out) :: A(nx,nx)
    integer :: i
    A        = 0.0_dp
    !A(1,1)   = 1.0_dp
    !A(nx,nx) = 1.0_dp
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
    !   --                                                                                                                                            --      --     --     --      --
    A(1,1)   =  1.0_dp  -  u * dt / dx - 2.0_dp * kappa * dt / (dx*dx)
    A(nx,nx) =  1.0_dp  -  u * dt / dx - 2.0_dp * kappa * dt / (dx*dx)

    do i = 2, nx-1
      A(i,i-1)   = kappa * dt / (dx*dx)        
      A(i,i  )   = 1.0_dp  -  u * dt / dx - 2.0_dp * kappa * dt / (dx*dx)
      A(i,i+1)   = u * dt / dx  +  kappa * dt / (dx*dx)

    end do


    A(1 ,nx  ) =  u * dt / dx  +  kappa * dt / (dx*dx)  
    A(1 ,1   ) =  1.0_dp  -  u * dt / dx - 2.0_dp * kappa * dt / (dx*dx)
    A(1 ,2   ) =  kappa * dt / (dx*dx)    

    A(nx,nx-1 ) = u * dt / dx  +  kappa * dt / (dx*dx) 
    A(nx,nx   ) = 1.0_dp  -  u * dt / dx - 2.0_dp * kappa * dt / (dx*dx)
    A(nx,1    ) = kappa * dt / (dx*dx)        


  end subroutine build_A

  !============================================================
  subroutine make_truth_ic(nx, dx, x0)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx
    real(dp), intent(out) :: x0(nx)
    integer :: i
    do i = 1, nx
      x0(i) = exp( - (( (real(i-1,dp)*dx) - 0.25_dp )**2) / (2.0_dp*(0.05_dp**2)) ) + &
              0.5_dp * exp( - (( (real(i-1,dp)*dx) - 0.65_dp )**2) / (2.0_dp*(0.08_dp**2)) )
    end do
  end subroutine make_truth_ic

  !============================================================
  function identity(n) result(Id)
    integer, intent(in) :: n
    real(dp) :: Id(n,n)
    integer :: ii
    Id = 0.0_dp
    do ii = 1, n
      Id(ii,ii) = REAL(1.0_dp)
    end do
  end function identity

  !============================================================
  function randn() result(r)
    real(dp) :: r
    real(dp) :: u1,u2
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0e-12_dp) u1 = 1.0e-12_dp
    r = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  !============================================================
  function randn_vec(n) result(v)
    integer, intent(in) :: n
    real(dp) :: v(n)
    integer :: i
    do i = 1, n
      v(i) = randn()
    end do
  end function randn_vec

  !============================================================
  subroutine cg_solve(A, b, x, maxit, tol)
    ! Conjugate Gradients for symmetric positive-definite A
    real(dp), intent(in) :: A(:,:), b(:)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: maxit
    real(dp), intent(in) :: tol
    integer :: n, i, it
    real(dp), allocatable :: r(:), p(:), Ap(:)
    real(dp) :: alpha, beta, rsold, rsnew, pAp

    n = size(b)
    allocate(r(n), p(n), Ap(n))
    r = b - matmul(A, x)
    p = r
    rsold = dot_product(r,r)
    do it = 1, maxit
      Ap = matmul(A, p)
      pAp = dot_product(p, Ap)
      if (abs(pAp) < 1.0e-30_dp) exit
      alpha = rsold / pAp
      x = x + alpha * p
      r = r - alpha * Ap
      rsnew = dot_product(r,r)
      if (sqrt(rsnew) < tol) exit
      beta = rsnew / rsold
      p = r + beta * p
      rsold = rsnew
    end do
    deallocate(r, p, Ap)
  end subroutine cg_solve

  !============================================================
  function rmse_vec(a, b) result(rmse)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: rmse
    integer :: n
    n = size(a)
    rmse = sqrt(sum((a - b)**2) / real(n, dp))
  end function rmse_vec

end program advec_diff_4dvar
