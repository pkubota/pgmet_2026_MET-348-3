! enkf_advec_diff.f90
program enkf_advec_diff
  implicit none
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: nx = 60, nt = 120
  integer, parameter :: Ne = 40            ! ensemble size
  real(dp), parameter :: L = 1.0_dp, u = 1.0_dp, kappa = 0.01_dp
  real(dp), parameter :: dx = L/(nx-1), dt = 0.25_dp*dx
  real(dp), parameter :: sigma_b = 0.10_dp, corr_len = 0.15_dp
  real(dp), parameter :: sigma_obs = 0.02_dp
  real(dp), parameter :: inflation = 1.02_dp   ! multiplicative inflation
  logical, parameter :: localize = .true.
  real(dp), parameter :: loc_radius = 0.25_dp ! localization radius (in domain units)

  integer, parameter :: m = 6
  integer, parameter :: obs_idx(m) = (/ 6,14,22,30,38,46 /)

  real(dp), allocatable :: A(:,:), truth(:,:), X(:,:,:), Xf(:,:,:), Xa_mean(:,:), y(:,:)
  real(dp), allocatable :: R(:,:), Rinv(:,:), H(:,:), Ht(:,:)
  real(dp), allocatable :: loc(:,:), x_mean(:)
  real(dp), allocatable :: Xp(:,:)
  real(dp), allocatable :: Yo(:,:)
  real(dp), allocatable :: Pf(:,:)
  real(dp), allocatable :: S(:,:), Sinv(:,:), K(:,:)

  integer :: i,j,kk,t, outunit, rmseunit

  call seed_rng(22222)

  allocate(A(nx,nx))
  allocate(truth(nt+1,nx))
  allocate(X(nx,Ne,nt+1))   ! X(:,j,t): ensemble member j at time t
  allocate(Xf(nx,Ne,nt+1))
  allocate(Xa_mean(nt+1,nx))
  allocate(y(nt,m))
  allocate(R(m,m), Rinv(m,m))
  allocate(H(m,nx), Ht(nx,m))
  allocate(loc(nx,nx))
  allocate(x_mean(nx))
  allocate(Yo(m,Ne))
  allocate(Xp(nx,Ne))
  allocate(Pf(nx,nx))
  allocate(S(m,m), Sinv(m,m), K(nx,m))

  call build_A(A, nx, dx, dt, u, kappa)

  ! truth
  call make_truth_ic(nx, dx, truth(1,1:nx))
  do t = 1, nt
    truth(t+1,1:nx) = matmul(A, truth(t,1:nx))
  end do

  ! observations
  y=0.0_dp
  do t = 1, nt
    do i = 1, m
      y(t,i) = truth(t+1, obs_idx(i)) + sigma_obs * randn()
    end do
  end do
  ! initial ensemble around truth(1)
  do j = 1, Ne
    do i = 1, nx
      X(i,j,1) = truth(1,i) + sigma_b * randn()
    end do
  end do
  Xa_mean(1,1:nx) = sum(X(1:nx,1:Ne,1), dim=2) / real(Ne,dp)

  ! R, H
  R = 0.0_dp; Rinv = 0.0_dp
  do i = 1, m
    R(i,i) = sigma_obs**2; Rinv(i,i) = 1.0_dp / R(i,i)
  end do
  H = 0.0_dp
  do i = 1, m
    H(i, obs_idx(i)) = 1.0_dp
  end do
  Ht = transpose(H)

  ! localization matrix (Gaspari-Cohn)
  if (localize) call build_localization(loc, nx, dx, loc_radius)

  open(newunit=outunit, file="enkf_out.dat", status="replace")
  open(newunit=rmseunit, file="enkf_rmse.dat", status="replace")
  write(outunit,'(A)') '# t truth(...) analysis_mean(...)'
  write(rmseunit,'(A)') '# t rmse_forecast rmse_analysis'

  ! EnKF cycle
  do t = 1, nt
    ! forecast ensemble: Xf(:,j) = A * X(:,j)
    do j = 1, Ne
      X(1:nx,j,t+1) = matmul(A, X(1:nx,j,t))
    end do

    ! compute forecast mean and perturbations
    x_mean(1:nx) = 0.0_dp
    do j = 1, Ne
      x_mean(1:nx) = x_mean(1:nx) + X(1:nx,j,t+1)
    end do
    x_mean(1:nx) = x_mean(1:nx) / real(Ne,dp)

    ! perturb observations matrix Y (m x Ne)
    do j = 1, Ne
      do i = 1, m
        Yo(i,j) = y(t,i) + sigma_obs * randn()
      end do
    end do

    ! build Pf ensemble (sample covariance) (nx x nx)
    do j = 1, Ne
      Xp(1:nx,j) = X(1:nx,j,t+1) - x_mean(1:nx)
    end do

    ! optionally apply inflation
    do j = 1, Ne
      X (1:nx,j,t+1) = x_mean(1:nx)  + inflation*(X(1:nx,j,t+1) - x_mean(1:nx))
      Xp(1:nx,j)     = X(1:nx,j,t+1) - x_mean(1:nx)
    end do

    ! compute sample covariance Pf = (1/(Ne-1)) Xp Xp^T
    Pf(1:nx,1:nx) = matmul(Xp(1:nx,1:Ne), transpose(Xp(1:nx,1:Ne))) / real(Ne-1,dp)

    ! apply localization if requested
    if (localize) Pf(1:nx,1:nx) = Pf(1:nx,1:nx) * loc(1:nx,1:nx)

    ! compute Kalman gain K = Pf H^T (H Pf H^T + R)^{-1}
    S(1:m,1:m) = matmul(H(1:m,1:nx), matmul(Pf(1:nx,1:nx), transpose(H(1:m,1:nx)))) + R(1:m,1:m)

    call invert_matrix(S(1:m,1:m), Sinv(1:m,1:m), m)

    K(1:nx,1:m) = matmul(Pf(1:nx,1:nx), matmul(transpose(H(1:m,1:nx)), Sinv(1:m,1:m)))

    ! update each ensemble member (stochastic EnKF)
    do j = 1, Ne
      ! ensemble member j forecast ,Ne
      X(1:nx,j,t+1) = X(1:nx,j,t+1) + matmul(K(1:nx,1:m), (Yo(1:m,j) - matmul(H(1:m,1:nx), X(1:nx,j,t+1))))
    end do

    ! analysis mean
    Xa_mean(t+1,1:nx) = sum(X(1:nx,1:m,t+1), dim=2) / real(Ne,dp)

    ! write time output
    write(outunit,'(I6,1X,60F12.6,1X,60F12.6)') t, truth(t+1,:), Xa_mean(t+1,1:nx)
    write(rmseunit,'(I6,2F12.6)') t, rmse_vec( sum(X(:,:,t+1), dim=2)/real(Ne,dp), truth(t+1,:) ), 0.0_dp

  end do
    deallocate(Yo, Xp, Pf, S, Sinv, K)

  close(outunit); close(rmseunit)
  print *, 'EnKF done. Outputs: enkf_out.dat, enkf_rmse.dat'

contains

  subroutine build_A(A, nx, dx, dt, u, kappa)
    integer, intent(in):: nx
    real(dp), intent(in):: dx, dt, u, kappa
    real(dp), intent(out):: A(nx,nx)
    integer :: i
    A = 0.0_dp; A(1,1) = 1.0_dp; A(nx,nx)=1.0_dp
    do i = 2, nx-1
      A(i,i-1) = u*dt/dx + kappa*dt/(dx*dx)
      A(i,i)   = 1.0_dp - u*dt/dx - 2.0_dp*kappa*dt/(dx*dx)
      A(i,i+1) = kappa*dt/(dx*dx)
    end do
  end subroutine build_A

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

  subroutine build_localization(rho, nx, dx, radius)
    real(dp), intent(out) :: rho(nx,nx)
    integer, intent(in) :: nx
    real(dp), intent(in) :: dx, radius
    integer :: i,j
    real(dp) :: xi, xj, d, r
    do i = 1, nx
      xi = real(i-1,dp)*dx
      do j = 1, nx
        xj = real(j-1,dp)*dx
        d = abs(xi - xj)
        r = d / radius
        rho(i,j) = gaspari_cohn(r)
      end do
    end do
  end subroutine build_localization

  function gaspari_cohn(r) result(w)
    real(dp), intent(in) :: r
    real(dp) :: w
    if (r >= 2.0_dp) then
      w = 0.0_dp
    else if (r <= 1.0_dp) then
      w = 1.0_dp - 5.0_dp*r**2/3.0_dp + 5.0_dp*r**3/4.0_dp + 0.5_dp*r**4 - 1.0_dp*r**5/4.0_dp
    else
      w = 4.0_dp - 5.0_dp*r + 5.0_dp*r**2/3.0_dp + 5.0_dp*r**3/8.0_dp - 0.5_dp*r**4 + 1.0_dp*r**5/12.0_dp
      w = w / r**2
    end if
  end function gaspari_cohn

  subroutine invert_matrix(A, Ainv, n)
    real(dp), intent(in) :: A(n,n)
    real(dp), intent(out) :: Ainv(n,n)
    integer, intent(in) :: n
    real(dp), allocatable :: aug(:,:)
    integer :: i,j,k,piv
    real(dp) :: maxv,fac,tmp
    allocate(aug(n,2*n))
    aug(:,1:n) = A; 
    aug(:,n+1:2*n) = 0.0_dp
    do i=1,n; 
       aug(i,n+i)=1.0_dp; 
    end do
    do i=1,n
      piv = i; 
      maxv = abs(aug(i,i))
      do k=i+1,n
        if (abs(aug(k,i)) > maxv) then; 
           maxv=abs(aug(k,i));
           piv=k; 
        end if
      end do
      if (piv /= i) then 
         aug([i,piv],:) = aug([piv,i],:); 
      end if
      if (abs(aug(i,i)) < 1.0e-16_dp) aug(i,i) = 1.0e-16_dp
      fac = aug(i,i); aug(i,:) = aug(i,:)/fac
      do j=1,n
        if (j /= i) then 
           tmp = aug(j,i); 
           aug(j,:) = aug(j,:) - tmp * aug(i,:); 
        end if
      end do
    end do
    Ainv = aug(:, n+1:2*n)
    deallocate(aug)
  end subroutine invert_matrix

  function randn() result(r)
    real(dp) :: r,u1,u2
    call random_number(u1); 
    call random_number(u2)
    if (u1 < 1.0e-16_dp) u1 = 1.0e-16_dp
    r = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  function rmse_vec(a,b) result(rmse)
    real(dp), intent(in) :: a(:), b(:)
    integer :: n; real(dp) :: rmse
    n = size(a); 
    rmse = sqrt(sum((a-b)**2) / real(n,dp))
  end function rmse_vec

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n,i; integer, allocatable :: s(:)
    call random_seed(size=n); allocate(s(n))
    do i=1,n;
       s(i)=i+seed; 
    end do
    call random_seed(put=s); 
    deallocate(s)
  end subroutine seed_rng

end program enkf_advec_diff
