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
! advec_diff_4dvar_precond.f90
! 4D-Var (strong-constraint) for 1D advection-diffusion with:
!  - Exponential B (non-diagonal)
!  - Preconditioned Conjugate Gradient (Jacobi preconditioner using diag(LHS))
!  - Configurable window length (win) and stride (stride)
!
program advec_diff_4dvar_precond
  use LinearSolve, only :inverse
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307)

  ! ---------------- User parameters ----------------
  integer, parameter :: nx = 60        ! number of grid points
  integer, parameter :: nt = 120       ! total time steps
  integer, parameter :: win = 12       ! window length (steps)
  integer, parameter :: stride = 6     ! stride between windows (stride <= win allowed)
  real(dp), parameter :: Ldom = 1.0_dp
  real(dp), parameter :: u = 1.0_dp
  real(dp), parameter :: kappa = 0.01_dp
  real(dp), parameter :: dx = Ldom / real(nx-1, dp)
  real(dp), parameter :: dt = 0.25_dp * dx   ! stability choice for explicit scheme

  real(dp), parameter :: sigma_b = 0.10_dp   ! background std
  real(dp), parameter :: corr_len = 0.15_dp  ! correlation length (for exponential B)
  real(dp), parameter :: sigma_obs = 0.02_dp ! observation std

  integer, parameter :: m = 6
  integer, parameter :: obs_idx(m) = (/ 6, 14, 22, 30, 38, 46 /)

  integer, parameter :: max_cg_iter = 2000
  real(dp), parameter :: cg_tol = 1.0e-8_dp

  ! ---------------- arrays ----------------
  real(dp), allocatable :: A(:,:), Ap(:,:,:), ATp(:,:,:)
  real(dp), allocatable :: truth(:,:), bkg(:,:), analysis(:,:)
  real(dp), allocatable :: x0_true(:), x0_b(:)
  real(dp), allocatable :: y(:,:)
  real(dp), allocatable :: B(:,:), Binv(:,:), R(:,:), Rinv(:,:), H(:,:), Ht(:,:)

  real(dp), allocatable :: LHS(:,:), LHS_inv(:,:), RHS(:), x0_opt(:)

  integer :: i,j,k,t,win_k,info
  integer :: outunit, rmseunit

  ! ---------------- init ----------------
  call seed_rng(13579)

  allocate(A(nx,nx))
  allocate(Ap(nx,nx,win+1))
  allocate(ATp(nx,nx,win+1))
  allocate(truth(nt+1,nx), bkg(nt+1,nx), analysis(nt+1,nx))
  allocate(x0_true(nx), x0_b(nx))
  allocate(y(nt,m))
  allocate(B(nx,nx), Binv(nx,nx), R(m,m), Rinv(m,m))
  allocate(H(m,nx), Ht(nx,m))
  allocate(LHS(nx,nx), LHS_inv(nx,nx), RHS(nx), x0_opt(nx))

  ! build model operator A
  call build_A(A, nx, dx, dt, u, kappa)

  ! compute powers A^k and (A^T)^k up to win
  Ap (:,:,1) = identity(nx)
  ATp(:,:,1) = identity(nx)

  do k = 2, win+1
    Ap (:,:,k)  = matmul(A, Ap(:,:,k-1))
    ATp(:,:,k) = transpose(Ap(:,:,k))
  end do

  ! truth IC and propagate truth

  call make_truth_ic(nx, dx, x0_true)

  truth(1,:) = x0_true

  do t = 1, nt
    truth(t+1,:) = matmul(A, truth(t,:))
  end do

  ! generate noisy observations at each time at obs_idx
  do t = 1, nt
    do i = 1, m
      y(t,i) = truth(t+1, obs_idx(i)) + sigma_obs * randn()
    end do
  end do

  ! background initial condition and free run
  x0_b = x0_true + sigma_b * randn_vec(nx)
  bkg(1,:) = x0_b
  do t = 1, nt
    bkg(t+1,:) = matmul(A, bkg(t,:))
  end do

  ! build exponential B and its inverse (dense)
  call build_B_exponential(B, nx, dx, sigma_b, corr_len)
  call invert_matrix(B, Binv, nx)

  ! build R and its inverse
  R = 0.0_dp
  do i = 1, m
    R(i,i) = sigma_obs**2
    Rinv(i,i) = 1.0_dp / R(i,i)
  end do

  ! build H (selection) and Ht
  H = 0.0_dp
  do i = 1, m
    H(i, obs_idx(i)) = 1.0_dp
  end do
  Ht = transpose(H)

  ! open outputs
  open(newunit=outunit, file="4dvar_out.dat", status="replace", action="write")
  write(outunit, '(A)') '# t truth(...) background(...) analysis(...)'
  open(newunit=rmseunit, file="rmse_4dvar.dat", status="replace", action="write")
  write(rmseunit, '(A)') '# t rmse_background rmse_analysis'

  ! initialize analysis with background
  analysis(1,:) = bkg(1,:)

  ! ---------- WINDOW LOOP WITH STRIDE ----------
  t = 1
  do while (t <= nt)
    ! determine actual window length available
    win_k = min(win, nt - t + 1)

    ! Build LHS = Binv + sum_{k=1..win_k} AT^k * H^T * Rinv * H * A^k
    LHS = Binv
    do k = 1, win_k
      LHS = LHS + matmul( ATp(:,:,k+1), matmul( matmul(Ht, Rinv), matmul(H, Ap(:,:,k+1)) ) )
    end do

    ! Build RHS = Binv * xb0 + sum_k AT^k * H^T * Rinv * y_k
    RHS = 0.0_dp
    RHS = matmul(Binv, bkg(t,:))
    do k = 1, win_k
      RHS = RHS + matmul( ATp(:,:,k+1), matmul(Ht, matmul(Rinv, y(t + k - 1, :))) )
    end do

    ! Solve LHS * x0_opt = RHS using direct matrix inversion
    call inverse(LHS, LHS_inv, nx)
    x0_opt = matmul(LHS_inv, RHS)

    ! propagate optimal x0 through window and write into analysis
    do k = 0, win_k-1
      analysis(t + k, :) = matmul(Ap(:,:,k+1), x0_opt)
    end do
    analysis(t + win_k, :) = matmul(Ap(:,:,win_k+1), x0_opt)

    ! advance by stride (allow overlap if stride < win)
    t = t + stride
  end do

  ! write outputs (truth, background, analysis) and RMSE
  do t = 1, nt+1
    write(outunit, '(I6,1X,60F12.6,1X,60F12.6,1X,60F12.6)') t-1, truth(t,:), bkg(t,:), analysis(t,:)
    write(rmseunit, '(I6,2F12.6)') t-1, rmse_vec(bkg(t,:), truth(t,:)), rmse_vec(analysis(t,:), truth(t,:))
  end do

  close(outunit)
  close(rmseunit)
  print *, 'Finished. Outputs: 4dvar_out.dat , rmse_4dvar.dat'

contains

  ! ---------------- utility: identity ----------------
  function identity(n) result(Id)
    integer, intent(in) :: n
    real(dp) :: Id(n,n)
    integer :: i
    Id = 0.0_dp
    do i = 1, n
      Id(i,i) = 1.0_dp
    end do
  end function identity

  ! ---------------- build A (explicit adv-diff) ----------------

  subroutine build_A(A, nx, dx, dt, u, kappa)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt, u, kappa
    real(dp), intent(out) :: A(nx,nx)
    integer :: i
    A        = 0.0_dp
    !A(1,1)   = 1.0_dp
    !A(nx,nx) = 1.0_dp
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

  ! ---------------- truth IC ----------------
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

  ! ---------------- build exponential B ----------------
  subroutine build_B_exponential(B, nx, dx, sigma_b, corr_len)
    real(dp), intent(out) :: B(nx,nx)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, sigma_b, corr_len
    integer :: i,j
    real(dp) :: xi, xj, dist
    do i = 1, nx
      xi = (real(i-1,dp) * dx)
      do j = 1, nx
        xj = (real(j-1,dp) * dx)
        dist = abs(xi - xj)
        B(i,j) = sigma_b**2 * exp(-dist / corr_len)
      end do
    end do
  end subroutine build_B_exponential

  ! ---------------- matrix inversion (Gauss-Jordan pivoting) ----------------
  subroutine invert_matrix(A, Ainv, n)
    real(dp), intent(in) :: A(n,n)
    real(dp), intent(out) :: Ainv(n,n)
    integer, intent(in) :: n
    real(dp), allocatable :: aug(:,:)
    integer :: i,j,k, pivrow
    real(dp) :: maxv, factor, tmp
    allocate(aug(n,2*n))
    aug(:,1:n) = A
    aug(:,n+1:2*n) = 0.0_dp
    do i = 1, n
      aug(i,n+i) = 1.0_dp
    end do

    do i = 1, n
      ! pivot selection (partial pivoting)
      pivrow = i
      maxv = abs(aug(i,i))
      do k = i+1, n
        if (abs(aug(k,i)) > maxv) then
          maxv = abs(aug(k,i)); pivrow = k
        end if
      end do
      if (pivrow /= i) then
        aug([i, pivrow], :) = aug([pivrow, i], :)
      end if
      if (abs(aug(i,i)) < 1.0e-16_dp) then
        print *, 'invert_matrix: near-zero pivot'
        aug(i,i) = 1.0e-16_dp
      end if
      ! normalize row i
      factor = aug(i,i)
      aug(i,:) = aug(i,:) / factor
      ! eliminate others
      do j = 1, n
        if (j /= i) then
          tmp = aug(j,i)
          aug(j,:) = aug(j,:) - tmp * aug(i,:)
        end if
      end do
    end do

    Ainv = aug(:, n+1:2*n)
    deallocate(aug)
  end subroutine invert_matrix

  ! ---------------- PCG with Jacobi preconditioner ----------------
  subroutine pcg_jacobi(A, b, x, Minv, maxit, tol)
    real(dp), intent(in) :: A(:,:), b(:), Minv(:)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: maxit
    real(dp), intent(in) :: tol
    integer :: n, it
    real(dp), allocatable :: r(:), z(:), p(:), Ap(:)
    real(dp) :: alpha, beta, rz_old, rz_new, pAp

    n = size(b)
    allocate(r(n), z(n), p(n), Ap(n))

    r = b - matmul(A, x)
    ! z = Minv * r (elementwise)
    z = r * Minv
    p = z
    rz_old = dot_product(r, z)

    do it = 1, maxit
      Ap = matmul(A, p)
      pAp = dot_product(p, Ap)
      if (abs(pAp) < 1.0e-30_dp) exit
      alpha = rz_old / pAp
      x = x + alpha * p
      r = r - alpha * Ap
      if (sqrt(dot_product(r,r)) < tol) exit
      z = r * Minv
      rz_new = dot_product(r, z)
      beta = rz_new / rz_old
      p = z + beta * p
      rz_old = rz_new
    end do

    deallocate(r, z, p, Ap)
  end subroutine pcg_jacobi

  ! ---------------- RNG (Box-Muller) ----------------
  function randn() result(r)
    real(dp) :: r
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    if (u1 < 1.0e-16_dp) u1 = 1.0e-16_dp
    r = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  function randn_vec(n) result(v)
    integer, intent(in) :: n
    real(dp) :: v(n)
    integer :: i
    do i = 1, n
      v(i) = randn()
    end do
  end function randn_vec

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: s(:)
    call random_seed(size=n)
    allocate(s(n))
    do i = 1, n
      s(i) = i + seed
    end do
    call random_seed(put=s)
    deallocate(s)
  end subroutine seed_rng

  ! ---------------- RMSE ----------------
  function rmse_vec(a, b) result(rmse)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: rmse
    integer :: n
    n = size(a)
    rmse = sqrt(sum((a - b)**2) / real(n, dp))
  end function rmse_vec

end program advec_diff_4dvar_precond
