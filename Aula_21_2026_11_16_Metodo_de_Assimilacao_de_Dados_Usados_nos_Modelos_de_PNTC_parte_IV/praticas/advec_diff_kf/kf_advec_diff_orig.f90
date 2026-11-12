! kf_advec_diff.f90
program kf_advec_diff
  implicit none
  integer , parameter :: dp = selected_real_kind(15,307)
  integer , parameter :: nx = 60, nt = 120
  real(dp), parameter :: L = 1.0_dp, u = 1.0_dp, kappa = 0.01_dp
  real(dp), parameter :: dx = L/(nx-1), dt = 0.25_dp*dx
  real(dp), parameter :: sigma_b = 0.10_dp, corr_len = 0.15_dp
  real(dp), parameter :: sigma_obs = 0.02_dp, q_diag = 1.0e-5_dp

  integer, parameter :: m = 6
  integer, parameter :: obs_idx(m) = (/ 6,14,22,30,38,46 /)

  real(dp), allocatable :: A(:,:), truth(:,:), forecast(:,:), analysis(:,:)
  real(dp), allocatable :: x_true(:), x_f(:), x_a(:)
  real(dp), allocatable :: P(:,:), Pf(:,:), Qa(:,:), Q(:,:), R(:,:), Rinv(:,:)
  real(dp), allocatable :: H(:,:), Ht(:,:)
  real(dp), allocatable :: B(:,:), Binv(:,:)
  real(dp), allocatable :: y(:,:)
  integer :: i,j,t, outunit, rmseunit

  call seed_rng(11111)

  allocate(A(nx,nx))
  allocate(truth(nt+1,nx), forecast(nt+1,nx), analysis(nt+1,nx))
  allocate(x_true(nx), x_f(nx), x_a(nx))
  allocate(P(nx,nx), Pf(nx,nx), Q(nx,nx), R(m,m), Rinv(m,m))
  allocate(H(m,nx), Ht(nx,m))
  allocate(B(nx,nx), Binv(nx,nx))
  allocate(y(nt,m))

  call build_A(A, nx, dx, dt, u, kappa)
  call build_B_exponential(B, nx, dx, sigma_b, corr_len)
  call invert_matrix(B, Binv, nx)

  ! build truth
  call make_truth_ic(nx, dx, x_true)
  truth(1,:) = x_true
  do t = 1, nt
    truth(t+1,:) = matmul(A, truth(t,:))
  end do

  ! observations
  do t = 1, nt
    do i = 1, m
      y(t,i) = truth(t+1, obs_idx(i)) + sigma_obs * randn()
    end do
  end do

  ! initial analysis x_a(0) = background (perturbed)
  x_a = x_true + sigma_b * randn_vec(nx)
  analysis(1,:) = x_a

  ! initial P = B
  P = B

  ! Q (model error) small diagonal
  Q = 0.0_dp
  do i = 1, nx
    Q(i,i) = q_diag
  end do

  ! H and R
  H = 0.0_dp
  do i = 1, m
    H(i, obs_idx(i)) = 1.0_dp
  end do
  Ht = transpose(H)
  R = 0.0_dp
  do i = 1, m
    R(i,i) = sigma_obs**2
    Rinv(i,i) = 1.0_dp / R(i,i)
  end do

  open(newunit=outunit, file="kf_out.dat", status="replace")
  open(newunit=rmseunit, file="kf_rmse.dat", status="replace")
  write(outunit,'(A)')  '# t truth(...) forecast(...) analysis(...)'
  write(rmseunit,'(A)') '# t rmse_forecast rmse_analysis'

  ! Loop time: forecast + update
  do t = 1, nt
    ! forecast
    x_f = matmul(A, x_a)
    Pf = matmul(A, matmul(P, transpose(A))) + Q
    forecast(t+1,:) = x_f

    ! Kalman gain
    ! S = H Pf H^T + R   (m x m)
    ! K = Pf H^T S^{-1}  (nx x m)
    call compute_K_and_update(Pf, H, Ht, R, Rinv, x_f, y(t,:), x_a, P)

    analysis(t+1,:) = x_a

  end do

  ! write outputs & rmse
  do t = 1, nt+1
    write(outunit,'(I6,1X,60F12.6,1X,60F12.6,1X,60F12.6)') t-1, truth(t,:), forecast(t,:), analysis(t,:)
    write(rmseunit,'(I6,2F12.6)') t-1, rmse_vec(forecast(t,:), truth(t,:)), rmse_vec(analysis(t,:), truth(t,:))
  end do

  close(outunit); close(rmseunit)
  print *, 'KF done. Outputs: kf_out.dat, kf_rmse.dat'

contains
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

  subroutine build_A(A, nx, dx, dt, u, kappa)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, dt, u, kappa
    real(dp), intent(out) :: A(nx,nx)
    integer :: i
    A = 0.0_dp
    A(1,1) = 1.0_dp; A(nx,nx) = 1.0_dp
    do i = 2, nx-1
      A(i,i-1) = u*dt/dx + kappa*dt/(dx*dx)
      A(i,i)   = 1.0_dp - u*dt/dx - 2.0_dp*kappa*dt/(dx*dx)
      A(i,i+1) = kappa*dt/(dx*dx)
    end do
  end subroutine build_A

  subroutine build_B_exponential(B, nx, dx, sigma_b, corr_len)
    real(dp), intent(out) :: B(nx,nx)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, sigma_b, corr_len
    integer :: i,j
    real(dp) :: xi, xj, d
    do i = 1, nx
      xi = (real(i-1,dp) * dx)
      do j = 1, nx
        xj = (real(j-1,dp) * dx)
        d = abs(xi-xj)
        B(i,j) = sigma_b**2 * exp(-d/corr_len)
      end do
    end do
  end subroutine build_B_exponential

  subroutine compute_K_and_update(Pf, H, Ht, R, Rinv, xf, y_obs, xa, Pa)
    real(dp), intent(in) :: Pf(nx,nx), H(:,:), Ht(:,:), R(:,:), Rinv(:,:)
    real(dp), intent(in) :: xf(nx), y_obs(:)
    real(dp), intent(inout) :: xa(nx), Pa(nx,nx)
    real(dp) :: S(size(R,1),size(R,2)), Sinv(size(R,1),size(R,2))
    real(dp) :: K(nx,size(R,1)), innov(size(R,1))
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
    real(dp) :: Id(n,n)
    integer :: ii
    Id = 0.0_dp
    do ii=1,n; Id(ii,ii)=1.0_dp; end do
  end function identity
  ! ---------------- matrix inversion (Gauss-Jordan pivoting) ----------------

  subroutine invert_matrix(A, Ainv, n)
    real(dp), intent(in) :: A(n,n)
    real(dp), intent(out) :: Ainv(n,n)
    integer, intent(in) :: n
    real(dp), allocatable :: aug(:,:)
    integer :: i,j,k,piv
    real(dp) :: maxv,fac,tmp
    allocate(aug(n,2*n))
    aug(:,1:n) = A
    aug(:,n+1:2*n) = 0.0_dp
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    do i=1,n
      ! pivot selection (partial pivoting)
      piv = i; maxv = abs(aug(i,i))
      do k=i+1,n
        if (abs(aug(k,i)) > maxv) then; maxv = abs(aug(k,i)); piv = k; end if
      end do
      if (piv /= i) then
        aug([i,piv],:) = aug([piv,i],:)
      end if
      if (abs(aug(i,i)) < 1.0e-16_dp) aug(i,i) = 1.0e-16_dp
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
    real(dp) :: r, u1,u2
    call random_number(u1); call random_number(u2)
    if (u1 < 1.0e-16_dp) u1 = 1.0e-16_dp
    r = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  function randn_vec(n) result(v)
    integer, intent(in):: n
    real(dp) :: v(n); integer ii
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
    real(dp), intent(in):: a(:), b(:); real(dp) :: rmse; integer n
    n = size(a); rmse = sqrt(sum((a-b)**2)/real(n,dp))
  end function rmse_vec

end program kf_advec_diff
