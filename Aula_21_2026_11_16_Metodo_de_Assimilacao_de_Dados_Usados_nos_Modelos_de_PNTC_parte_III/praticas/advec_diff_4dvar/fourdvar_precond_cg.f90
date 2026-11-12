program fourdvar_precond_cg
  implicit none
  integer , parameter :: dp = kind(1.0d0)
  integer , parameter :: nx=50, nt=50, nobs=5
  real(dp), parameter :: dx=1.0_dp/nx, dt=0.001_dp
  real(dp), parameter :: u=1.0_dp, kappa=0.01_dp
  real(dp), parameter :: sig_b=1.0_dp, Lc=0.1_dp
  integer             :: i, n, iter, max_iter

  real(dp), dimension(nx)      :: xb, x0, grad_v, v
  real(dp), dimension(nx,nx)   :: B, L
  real(dp), dimension(nx,nt)   :: x_fwd, lam
  real(dp), dimension(nobs,nt) :: yobs
  real(dp), dimension(nx)      :: p, r, Ap

  integer, dimension(nobs) :: obs_idx
  real(dp), parameter      :: Rinv = 1.0_dp/0.25_dp
  real(dp)                 :: alpha, beta, rr_old, rr_new

  !--- inicialização
  obs_idx = (/5,10,20,30,45/)
  xb      = 10.0_dp
  v       = 0.0_dp

  call build_B(B, sig_b, Lc)
  
  call cholesky(B, L)

  call build_observations(yobs, obs_idx)

  max_iter = 40

  ! ---- CG iterations in preconditioned space (v)
  ! r = grad(v), p = -r

  call compute_gradient(v, xb, L, x_fwd, lam, grad_v, yobs, obs_idx)

  r      =  grad_v
  p      = -r
  rr_old = dot_product(r,r)

  do iter=1,max_iter
      call apply_Hessian_action(v, xb, L, p, Ap, yobs, obs_idx)
      alpha  = rr_old / dot_product(p,Ap)
      v      = v + alpha*p
      r      = r + alpha*Ap
      rr_new = dot_product(r,r)
      beta   = rr_new/rr_old
      p      = -r + beta*p
      rr_old = rr_new
      if (mod(iter,5)==0) write(*,'(A,I3,A,F12.6)') "Iter ",iter," |res|^2 =", rr_new
  end do

  ! reconstruir x0
  x0 = xb + matmul(L,v)

  write(*,*) "4D-Var concluido. Estado inicial otimo:"
  do i=1,nx
     write(*,'(F10.5)') x0(i)
  end do


contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine forward_model(x_init, x)
    real(dp), intent(in) :: x_init(nx)
    real(dp), intent(out):: x(nx,nt)
    real(dp):: c, d
    integer :: i,n
    x(:,1)=x_init
    c=u*dt/dx
    d=kappa*dt/(dx*dx)
    do n=1,nt-1
      do i=2,nx-1
        x(i,n+1)=x(i,n) -c*(x(i,n)-x(i-1,n))+d*(x(i+1,n)-2*x(i,n)+x(i-1,n))
      end do
      x(1 ,n+1)=x(2   ,n+1)
      x(nx,n+1)=x(nx-1,n+1)
    end do
  end subroutine forward_model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine adjoint_model(x, yobs, obs_idx, lam)
    real(dp), intent(in ) :: x(nx,nt), yobs(nobs,nt)
    integer , intent(in ) :: obs_idx(nobs)
    real(dp), intent(out) :: lam(nx,nt)
    real(dp):: c,d
    integer :: i,n,j

    lam=0.0_dp

    do n=1,nt
      do j=1,nobs
        lam(obs_idx(j),n) = lam(obs_idx(j),n) + Rinv*(yobs(j,n) - x(obs_idx(j),n))
      end do
    end do

    c=u*dt/dx
    d=kappa*dt/(dx*dx)

    do n=nt-1,1,-1
      do i=2,nx-1
        lam(i,n)=lam(i,n) + c*(lam(i+1,n+1)-lam(i,n+1)) &
                          + d*(lam(i+1,n+1)-2*lam(i,n+1)+lam(i-1,n+1))
      end do
      lam(1,n)=lam(2,n)
      lam(nx,n)=lam(nx-1,n)
    end do

  end subroutine adjoint_model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine build_B(B, sig, Lc)
    real(dp), intent(out):: B(nx,nx)
    real(dp), intent(in) :: sig,Lc
    integer :: i,j
    do i=1,nx
      do j=1,nx
        B(i,j)=sig*exp(-abs((i-j)*dx)/Lc)
      end do
    end do
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cholesky(A,L)
    real(dp), intent(in) :: A(nx,nx)
    real(dp), intent(out):: L(nx,nx)
    integer:: i,j,k
    L=0.0_dp
    do i=1,nx
      do j=1,i
        L(i,j)=A(i,j)
        do k=1,j-1
          L(i,j)=L(i,j)-L(i,k)*L(j,k)
        end do
        if(i==j) then
          L(i,j)=sqrt(L(i,j))
        else
          L(i,j)=L(i,j)/L(j,j)
        end if
      end do
    end do
  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine build_observations(yobs, obs_idx)
    real(dp), intent(out):: yobs(nobs,nt)
    integer, intent(in) :: obs_idx(nobs)
    integer :: i,n
    do n=1,nt
      do i=1,nobs
        yobs(i,n)=12.0_dp + 0.5_dp*sin(0.1_dp*n) + 0.1_dp*i
      end do
    end do
  end subroutine build_observations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_gradient(v,xb,L,x_fwd,lam,grad_v,yobs,obs_idx)
    real(dp), intent(in):: v(nx), xb(nx), L(nx,nx)
    real(dp), intent(out):: grad_v(nx)
    real(dp), intent(out):: x_fwd(nx,nt), lam(nx,nt)
    real(dp), intent(in):: yobs(nobs,nt)
    integer, intent(in):: obs_idx(nobs)

    call forward_model(xb+matmul(L,v), x_fwd)

    call adjoint_model(x_fwd, yobs   , obs_idx, lam)

    grad_v = v - matmul(L, lam(:,1))

  end subroutine compute_gradient
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine apply_Hessian_action(v, xb, L, p, Hp, yobs, obs_idx)
    real(dp), intent(in) :: v(nx), xb(nx), L(nx,nx), p(nx)
    real(dp), intent(out):: Hp(nx)
    real(dp), intent(in) :: yobs(nobs,nt)
    integer , intent(in) :: obs_idx(nobs)
    real(dp), dimension(nx) :: vp
    real(dp), dimension(nx,nt) :: xf, lam_temp
    vp = matmul(L,p)
    call forward_model(xb+matmul(L,v), xf)
    call adjoint_model(xf            , yobs, obs_idx, lam_temp)
    Hp = p - matmul(L, lam_temp(:,1))
  end subroutine

end program fourdvar_precond_cg
