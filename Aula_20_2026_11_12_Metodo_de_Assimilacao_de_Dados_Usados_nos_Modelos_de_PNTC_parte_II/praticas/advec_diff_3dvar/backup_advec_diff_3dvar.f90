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
  use LinearSolve, only :inverse
  implicit none
  integer, parameter :: nx = 50, nt = 50
  real(8), parameter :: L = 1.0d0, dt = 0.001d0, u = 1.0d0, K = 0.01d0
  real(8) :: dx
  real(8), dimension(nx) :: T, Tnew, Tb, Ta, y, diff
  real(8), dimension(nx,nx) :: B, Ht, H, R, Id, tmp,mat_A,Mat_B, gain
  integer :: i, j,n,idim,jdim, nobs(1:2)

  dx   = L / (nx-1)
  !Ponto observado:
  nobs(1) = nx/2
  nobs(2) = nx/2-1
  ! --- Condicao inicial (gaussiana) ---
  do i = 1, nx
     T(i) = exp(-((i*dx - 0.5d0)**2)/0.01d0)
  end do
  Tb(1:nx) = T(1:nx)

  ! --- Matriz de covariancia B (identidade simplificada) ---
  B(1:nx,1:nx) = 0.0d0
  do i = 1, nx
     B(i,i) = 0.05d0
  end do

  ! --- Matriz de erro de observacao R ---
  R(1:nx,1:nx) = 0.0d0
  DO i=1,size(nobs,dim=1)
     R(nobs(i),nobs(i)) = 0.02d0
  END DO
  ! --- Operador de observacao H (observa apenas 1 ponto central) ---
  H (1:nx,1:nx)        = 0.0d0
  DO i=1,size(nobs,dim=1)
     H(nobs(i),nobs(i)) = 1.0d0
  END DO
  Ht(1:nx,1:nx)       = transpose(H)

  ! --- Simulação do modelo (background) ---
  do n = 1, nt
      call advec_diff_step(T, Tnew, u, K, dx, dt, nx)
      T = Tnew
  end do
  Tb = T   ! estado background
  ! --- Observação sintetica ---
  y = 0.0d0
  DO i=1,size(nobs,dim=1)
     y(nobs(i)) = Tb( nobs(i)  ) + 0.1d0   ! adiciona erro positivo
  END DO

  ! --- Aplicação do 3DVar ---
  tmp = matmul(H, matmul(B, Ht)) + R
  
  ! --- Inverse Operator tmp ---
  
  jdim=size(nobs,dim=1)
  idim=size(nobs,dim=1)
  DO j=1,jdim
     DO i=1,idim
         mat_A(i,j) =tmp (nobs(i),nobs(j))
         Mat_B(i,j) =tmp (nobs(i),nobs(j))
     END DO
  END DO

  call inverse(mat_A(1:idim,1:jdim),Mat_B(1:idim,1:jdim),idim)

  DO j=1,jdim
     DO i=1,idim
        Id (nobs(i),nobs(j))= mat_B(i,j) 
     END DO
  END DO

  gain = matmul(B, matmul(Ht, Id))   ! ganho de Kalman

  Ta = Tb + matmul(gain, (y - matmul(H, Tb)))
  ! --- Resultado ---
  DO i=1,size(nobs,dim=1)
    print *, "----------------"
    print *, "Ponto observado:", nobs(i)
    print *, "Background:     ", Tb(nobs(i))
    print *, "Observacao:     ", y(nobs(i))
    print *, "Analise:        ", Ta(nobs(i))
  END DO
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
       Tnew(i) = T(i) - (u*dt/(2.d0*dx))*(T(i+1)-T(i-1)) + (K*dt/(dx*dx))*(T(i+1)-2*T(i)+T(i-1))
    end do
    ! Condicaes de contorno periodicas
    Tnew(1) = T(nx-1)
    Tnew(nx) = T(2)
  end subroutine advec_diff_step

  function inv1x1(a) result(ainv)
    real(8), dimension(1,1), intent(in) :: a
    real(8), dimension(1,1) :: ainv
    ainv(1,1) = 1.d0 / a(1,1)
  end function inv1x1

end program advec_diff_3dvar
