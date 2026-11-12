program fourdvar_advecdiff
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: nx = 50, nt = 50, nobs = 5
  real(dp), parameter :: dx = 1.0_dp/nx, dt = 0.001_dp
  real(dp), parameter :: u = 1.0_dp, kappa = 0.01_dp
  integer :: i, n, k, iter, max_iter
  real(dp) :: alpha     ! passo do gradiente
  real(dp) :: c, d, Jcost

  ! Vetores de estado
  real(dp), dimension(nx)     :: x0, xb, gradJ
  real(dp), dimension(nx, nt) :: x_fwd, lambda

  ! Observações
  integer , dimension(nobs)      :: obs_idx
  real(dp), dimension(nobs, nt) :: yobs
  real(dp), dimension(nx)       :: tmp

  ! Covariâncias (aqui simplificadas como diagonais)
  real(dp), parameter :: Binv = 1.0_dp/1.0_dp
  real(dp), parameter :: Rinv = 1.0_dp/0.25_dp

  ! Inicialização

  xb(1:nx)      = 10.0_dp
  x0(1:nx)      = xb(1:nx)                     ! chute inicial da análise
  obs_idx = (/ 5, 10, 20, 30, 45 /)

  ! Criar campo com verdade + ruído como observação 
  do n = 1, nt
     do i = 1, nobs
        yobs(i,n) = 12.0_dp + 0.5_dp*sin(0.1_dp*n) + 0.1_dp*real(i,dp)
     end do
  end do

  max_iter = 10
  alpha    = 0.5_dp

  do iter = 1, max_iter

      ! 1) Modelo direto
      call forward_model(x0(1:nx)         , x_fwd(1:nx, 1:nt)) 

      ! 2) Modelo adjunto ? calcula lambda(0)
      call adjoint_model(x_fwd(1:nx, 1:nt), yobs(1:nobs, 1:nt), obs_idx, lambda(1:nx, 1:nt))

      ! 3) Gradiente do custo
      gradJ(1:nx) = Binv*(x0(1:nx) - xb(1:nx)) - lambda(1:nx,1)

      ! 4) Atualização (gradiente descendente)
      x0(1:nx) = x0(1:nx) - alpha*gradJ(1:nx)

      ! 5) Valor do custo (opcional)
      Jcost = 0.5_dp*sum((x0(1:nx)-xb(1:nx))**2)*Binv

      do n=1,nt
         tmp(1:nx) = x_fwd(1:nx,n)
         do i=1,nobs
            Jcost = Jcost + 0.5_dp*Rinv*( yobs(i,n) - tmp(obs_idx(i)) )**2
         end do
      end do 

      if (mod(iter,5)==0) write(*,'(A,I3,A,F15.6)') 'Iter=',iter,'  J=',Jcost

  end do

  write(*,*) '4D-Var concluido. Estado inicial otimo:'
  do i=1,nx
     write(*,'(F15.5)') x0(i)
  end do

contains

  subroutine forward_model(Xb_init, Xb)
    real(dp), intent(in) :: Xb_init(nx)
    real(dp), intent(out) :: Xb(nx,nt)
    real(dp) :: c, d
    integer :: i,n
    Xb(:,1) = Xb_init
    c = u*dt/dx
    d = kappa*dt/(dx*dx)
    do n=1,nt-1
       do i=2,nx-1
         Xb(i,n+1) = Xb(i,n) - c*(Xb(i,n)  -Xb(i-1,n)) &
                             + d*(Xb(i+1,n)-2*Xb(i,n)+Xb(i-1,n))
       end do
       Xb(1,n+1) = Xb(2,n+1)
       Xb(nx,n+1) = Xb(nx-1,n+1)
    end do
  end subroutine

  subroutine adjoint_model(x, yobs, obs_idx, lam)
    real(dp), intent(in) :: x(nx,nt), yobs(nobs,nt)
    integer , intent(in) :: obs_idx(nobs)
    real(dp), intent(out):: lam(nx,nt)
    real(dp) :: c, d
    integer :: n,i,j
    lam = 0.0_dp
    ! termo da derivada do mísmatch observacional
    do n = 1, nt
       do j = 1, nobs
          lam(obs_idx(j),n) = lam(obs_idx(j),n) + Rinv*(yobs(j,n)-x(obs_idx(j),n))
       end do
    end do
    ! integração backward
    c = u*dt/dx
    d = kappa*dt/(dx*dx)
    do n = nt-1,1,-1
       do i=2,nx-1
          lam(i,n) = lam(i,n) + c*(lam(i+1,n+1)-lam(i,n+1)) &
                              + d*(lam(i+1,n+1)-2*lam(i,n+1)+lam(i-1,n+1))
       end do
       lam(1,n) = lam(2,n)
       lam(nx,n) = lam(nx-1,n)
    end do
  end subroutine

end program fourdvar_advecdiff
