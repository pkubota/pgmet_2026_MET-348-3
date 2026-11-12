MODULE Class_WritetoGrads
 !USE Class_Fields, Only: Idim,xa
 IMPLICIT NONE
 PRIVATE
 INTEGER, PUBLIC      , PARAMETER  :: r8=8
 INTEGER, PUBLIC      , PARAMETER  :: r4=4
 INTEGER                    , PARAMETER :: UnitData=1
 INTEGER                    , PARAMETER :: UnitCtl=2
 CHARACTER (LEN=400)                   :: FileName
 LOGICAL                                            :: CtrlWriteDataFile
 INTEGER :: Idim
 REAL (kind=r4), ALLOCATABLE :: xa(:)
 PUBLIC :: SchemeWriteCtl
 PUBLIC :: SchemeWriteData
 PUBLIC :: InitClass_WritetoGrads
CONTAINS
 SUBROUTINE InitClass_WritetoGrads(Idim_in,xa_in)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Idim_in
    REAL(kind=r8) , INTENT(IN)   :: xa_in(:)
    ALLOCATE(xa(size(xa_in,dim=1)))
    Idim=Idim_in
    xa=xa_in
    FileName=''
    FileName='AdvecLinearConceitual1D'
    CtrlWriteDataFile=.TRUE.
 END SUBROUTINE InitClass_WritetoGrads

 FUNCTION SchemeWriteData(vars,irec)  RESULT (ok)
    IMPLICIT NONE
    REAL (KIND=r8), INTENT (INOUT) :: vars(Idim)
    INTEGER       , INTENT (INOUT) :: irec
    INTEGER        :: ok
    INTEGER        :: lrec
    REAL (KIND=r4) :: Yout(Idim)
    IF(CtrlWriteDataFile)INQUIRE (IOLENGTH=lrec) Yout
    IF(CtrlWriteDataFile)OPEN(UnitData,FILE=TRIM(FileName)//'.bin',&
    FORM='UNFORMATTED', ACCESS='DIRECT', STATUS='UNKNOWN', &
    ACTION='WRITE',RECL=lrec)
    ok=1
    CtrlWriteDataFile=.FALSE.
    Yout=REAL(vars(1:Idim),KIND=r4)
    irec=irec+1
    WRITE(UnitData,rec=irec)Yout
     ok=0
 END FUNCTION SchemeWriteData

 FUNCTION SchemeWriteCtl(nrec)  RESULT (ok)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nrec
    INTEGER             :: ok,i
    ok=1
   OPEN(UnitCtl,FILE=TRIM(FileName)//'.ctl',FORM='FORMATTED', &
   ACCESS='SEQUENTIAL',STATUS='UNKNOWN',ACTION='WRITE')
    WRITE (UnitCtl,'(A6,A           )')'dset ^',TRIM(FileName)//'.bin'
    WRITE (UnitCtl,'(A                 )')'title  EDO'
    WRITE (UnitCtl,'(A                 )')'undef  -9999.9'
    WRITE (UnitCtl,'(A6,I8,A8   )')'xdef  ',Idim,' levels '
    WRITE (UnitCtl,'(10F16.10   )')(xa(i),i=1,Idim)
    WRITE (UnitCtl,'(A                  )')'ydef  1 linear  -1.27 1'
    WRITE (UnitCtl,'(A6,I6,A25   )')'tdef  ',nrec,' linear  00z01jan0001 1hr'
    WRITE (UnitCtl,'(A20             )')'zdef  1 levels 1000 '
    WRITE (UnitCtl,'(A           )')'vars 3'
    WRITE (UnitCtl,'(A           )')'uc 0 99 resultado da edol yc'
    WRITE (UnitCtl,'(A           )')'ua 0 99 solucao analitica ya'
    WRITE (UnitCtl,'(A           )')'an 0 99 solucao analise   an'
    WRITE (UnitCtl,'(A           )')'endvars'
    CLOSE (UnitCtl,STATUS='KEEP')
    CLOSE (UnitData,STATUS='KEEP')
    ok=0
 END FUNCTION SchemeWriteCtl
END MODULE Class_WritetoGrads
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program advec_diff_4dvar_equiv
   use Class_WritetoGrads, only:SchemeWriteCtl,SchemeWriteData,InitClass_WritetoGrads
  !! 4D-Var forte para equacao 1D de adveccao-difusao
  !! Problema analogo a um KF 1D: mesmo tipo de verdade, observacoes, B, R e H.
  !!
  !! x_{n+1} = A x_n   (A: operador de um passo do modelo)
  !!
  !! J(x0) = 1/2 (x0 - xb)^T B^{-1} (x0 - xb)
  !!       + 1/2 sum_{k=1..win} [ y_k - H A^k x0 ]^T R^{-1} [ y_k - H A^k x0 ]
  !!
  !! Solucao: (linear forte) resolvendo LHS x0 = RHS, com
  !!   LHS = B^{-1} + sum_k (A^k)^T H^T R^{-1} H A^k
  !!   RHS = B^{-1} xb + sum_k (A^k)^T H^T R^{-1} y_k
  !!
  use iso_fortran_env, only: real64
  implicit none

  ! ------------------------
  ! parâmetros do problema
  ! ------------------------
  integer, parameter :: nx   = 60      ! pontos de grade
  integer, parameter :: nt   = 100      ! tempo total (numero de passos)
  integer, parameter :: win  = 10      ! janela de assimilacao (passos)
  integer, parameter :: nobs = 60      ! numero de pontos observados

  real(real64), parameter :: L     = 1.0_real64   ! comprimento do dominio [0,L]
  real(real64), parameter :: c_adv = 1.0_real64   ! velocidade de adveccao
  real(real64), parameter :: k_dif = 0.01_real64  ! coeficiente de difusao
  real(real64), parameter :: dt    = 0.001_real64  ! passo de tempo
  real(real64), parameter :: sigma_b   = 0.2_real64   ! desvio padrao do erro de background
  real(real64), parameter :: Lcorr     = 0.1_real64   ! comprimento de correlacao de B
  real(real64), parameter :: sigma_obs = 0.02_real64  ! desvio padrao do erro de observacao
  real(real64), parameter :: pi = 4.0_real64*atan(1.0_real64)

  ! ------------------------
  ! malha e estados
  ! ------------------------
  real(real64) :: dx
  real(real64), dimension(nx) :: x      ! coordenadas da malha
  real(real64), dimension(nx) :: x_true0, x_true      ! estado verdadeiro (inicial e instantaneo)
  real(real64), dimension(nx,0:nt) :: truth_traj      ! trajetaia verdadeira
  real(real64), dimension(nx,0:1) :: xf_bkg      ! trajetaia verdadeira
  real(real64), dimension(nx) :: x0_bkg, x_bkg        ! estado de background (inicial e instantaneo)
  real(real64), dimension(nx,0:nt) :: bkg_traj        ! trajetoria de background
  real(real64), dimension(nx,0:nt) :: x0t_ana               ! estado analisado inicial e instantaneo
  real(real64), dimension(nx) :: x0_ana, x_ana        ! estado analisado inicial e instantaneo
  real(real64), dimension(nx,0:nt) :: ana_traj        ! trajetoria analisada calculada pelo modelo

  ! ------------------------
  ! operador de modelo
  ! ------------------------
  real(real64), dimension(nx,nx)       :: A     ! operador de um passo (advec-dif)
  real(real64), dimension(nx,nx,0:win) :: Ap    ! A^k, k=0..win

  ! ------------------------
  ! matrizes de erro e operador de observação
  ! ------------------------
  real(real64), dimension(nx,nx)     :: B, Binv
  real(real64), dimension(nobs,nobs) :: R, Rinv
  real(real64), dimension(nobs,nx)   :: H
  real(real64), dimension(nx,nobs)   :: Ht

  ! ------------------------
  ! 4D-Var: LHS, RHS, etc.
  ! ------------------------
  real(real64), dimension(nx,nx) :: LHS, LHS_inv
  real(real64), dimension(nx)    :: RHS

  ! observacoes
  integer     , dimension(nobs)      :: obs_idx
  real(real64), dimension(nobs,0:nt) :: y_obs  ! y_obs(:,k) observacoes no tempo k

  ! utilitarios
  integer :: i, j, k, n, it, kk,irec,istat,icount
  real(real64) :: rand1, rand2, z, d
  real(real64) :: rmse_bkg, rmse_ana
  print*,'passei' 
  ! ------------------------
  ! inicializacao da malha
  ! ------------------------
  dx = L / real(nx,real64)
  do i = 1, nx
    x(i) = (real(i-1,real64) + 0.5_real64)*dx
  end do
  call InitClass_WritetoGrads(nx,x)
  ! ------------------------
  ! condição inicial verdadeira: por exemplo, sino gaussiano
  ! ------------------------
  call set_true_initial_condition(x, x_true0)

  ! construindo trajetoria verdadeira integrando o modelo
  call build_model_matrix(A, nx, dx, dt, c_adv, k_dif)

  call compute_powers_of_A(A, Ap, nx, win)

  call integrate_truth(A, x_true0, truth_traj, nx, nt)

  ! ------------------------
  ! background inicial = verdade + erro de background
  ! ------------------------
  call init_random_seed()
  call build_B_gaussian(B, Binv, x, nx, sigma_b, Lcorr)
  call sample_background(x_true0, x0_bkg, B, nx)

  ! trajetória de background propagando x0_bkg
  call integrate_truth(A, x0_bkg, bkg_traj, nx, nt)

  ! ------------------------
  ! definir pontos de observação (ex.: nobs pontos igualmente espacados)
  ! ------------------------
  call define_obs_points(obs_idx, nobs, nx)

  ! operador de observação H (seleção direta) e R
  call build_H_and_R(H, Ht, R, Rinv, obs_idx, nobs, nx, sigma_obs)

  ! gerar observações com ruído gaussiano
  DO k=0,nt
    call generate_observations(truth_traj(1:nx,k), y_obs(1:nx,k), H(1:nobs,1:nx), obs_idx, nobs, nx, sigma_obs)
  END DO
  ! ------------------------
  ! 4D-Var em uma janela [0..win]
  ! ------------------------

  x0t_ana(1:nx,0) = x_true0(1:nx)
  x0_bkg(1:nx)=x_true0(1:nx)
  n = win   ! tamanho da janela em número de passos

  do it = 1, nt

      ! LHS = B^{-1} + sum_{k=1..n} (A^k)^T H^T R^{-1} H A^k
      LHS(1:nx,1:nx) = Binv(1:nx,1:nx)
      kk=1 
      do k = MAX(it-int(win/2.0),1)+1, MIN(it+int(win/2.0)-1,nt)

         LHS(1:nx,1:nx) = LHS(1:nx,1:nx) + matmul( transpose(Ap(1:nx,1:nx,kk)), &
                                           matmul( matmul(Ht(1:nx,1:nobs), Rinv(1:nobs,1:nobs)), &
                                           matmul( H(1:nobs,1:nx), Ap(1:nx,1:nx,kk)) ) )
         kk=kk+1
      end do

      call integrate_truth(A, x0_bkg, xf_bkg(1:nx,0:1), nx, 1)
      x0_bkg(1:nx)=xf_bkg(1:nx,1)


      ! RHS = B^{-1} x0_bkg + sum_{k=1..n} (A^k)^T H^T R^{-1} y_k
      RHS(1:nx) = matmul(Binv(1:nx,1:nx), x0_bkg(1:nx))
      kk=1 
      do k = MAX(it-int(win/2.0),1)+1, MIN(it+int(win/2.0)-1,nt)
         ! gerar observações com ruído gaussiano
         call generate_observations(truth_traj(1:nx,k), y_obs(1:nobs,kk), H(1:nobs,1:nx), obs_idx(1:nobs), nobs, nx, sigma_obs)

         RHS(1:nx) = RHS(1:nx) + matmul( transpose(Ap(1:nx,1:nx,kk)), matmul(Ht(1:nx,1:nobs), &
                                 matmul(Rinv(1:nobs,1:nobs), y_obs(1:nobs,kk))) )
         kk=kk+1
      end do
      ! Resolver LHS x0_ana = RHS por inversão de matriz
      call invert_matrix(LHS(1:nx,1:nx), LHS_inv(1:nx,1:nx), nx)

      x0_ana(1:nx    ) = matmul(LHS_inv(1:nx,1:nx), RHS(1:nx))
      x0t_ana(1:nx,it) = x0_ana(1:nx)
      x0_bkg(1:nx    ) = x0_ana(1:nx)
  END DO

  ! ------------------------
  ! Avaliar RMSE background vs analise na janela
  ! ------------------------
  irec=0 
  icount=0
  istat=SchemeWriteData(bkg_traj  (1:nx,0),irec)
  istat=SchemeWriteData(truth_traj(1:nx,0),irec)
  istat=SchemeWriteData(ana_traj  (1:nx,0),irec)

  print *, 'Tempo   RMSE(bkg)   RMSE(ana)'
  do k = 1, nt
      ! Propagar condição inicial analisada
      call integrate_truth(A, x0t_ana(1:nx,k-1), ana_traj(:,0:1), nx, 1)

       do i = 1, nx
        print*,truth_traj(i,k),bkg_traj(i,k),ana_traj(i,1)
       end do
      istat=SchemeWriteData(bkg_traj  (1:nx,k),irec)
      istat=SchemeWriteData(truth_traj(1:nx,k),irec)
      istat=SchemeWriteData(ana_traj  (1:nx,1),irec)
      icount=icount+1
      print*, '==========================='
      call compute_rmse(truth_traj(:,k), bkg_traj(:,k), nx, rmse_bkg)
      call compute_rmse(truth_traj(:,k), ana_traj(:,1), nx, rmse_ana)
       !print '(I4,2F12.6)', k, rmse_bkg, rmse_ana
  end do
 istat=SchemeWriteCtl(icount)

contains

  !===========================
  subroutine set_true_initial_condition(x, x0)
    real(real64), intent(in)  :: x(:)
    real(real64), intent(out) :: x0(:)
    integer :: i, n
    real(real64) :: x0_center, width

    n = size(x)
    x0_center = 0.3_real64
    width     = 0.05_real64

    do i = 1, n
      x0(i) = exp( - ( (x(i)-x0_center)**2 ) / (2.0_real64*width**2) )
    end do
  end subroutine set_true_initial_condition
  !===========================

  !===========================
  subroutine build_model_matrix(A, nx, dx, dt, c, kappa)
    !! Constroi o operador de um passo A da equacao
    !!   du/dt + c du/dx = kappa d2u/dx2
    !!
    !! Esquema explicito simples com condicoes periodicas:
    !!   upwind para advecção + 2ª derivada centrada para difusao
    !!
    integer, intent(in)       :: nx
    real(real64), intent(out) :: A(nx,nx)
    real(real64) :: B(nx,nx)
    real(real64), intent(in)  :: dx, dt, c, kappa

    integer :: i, j
    real(real64) :: alpha, beta
    B = 0.0_real64

    A = 0.0_real64
    !
    ! A(1,1)   = 1.0_dp
    !
    ! A(nx,nx) = 1.0_dp
    !
    ! dA/dt + u dA/dx = kappa d2A/dx2
    ! dA/dt + u dA/dx = kappa d2A/dx2
    !
    ! A(n+1,i) = A(n,i) - (udt/dx)*A(n,i+1) + (udt/dx)*A(n,i)  + (k*dt/dxdx)*A(n,i+1) - (k*dt/dxdx)*2A(n,i) + (k*dt/dxdx)*A(n,i-1)
    !  
    ! A(n+1,i) = A(n,i) - (udt/dx)*A(n,i+1) + (udt/dx)*A(n,i)  + (k*dt/dxdx)*A(n,i+1) - (k*dt/dxdx)*2A(n,i) + (k*dt/dxdx)*A(n,i-1)
    !
    ! A(n+1,i) = (k*dt/dxdx)*A(n,i-1) + A(n,i) + (udt/dx)*A(n,i) - (k*dt/dxdx)*2A(n,i) - (udt/dx)*A(n,i+1)  + (k*dt/dxdx)*A(n,i+1)
    !
    ! A(n+1,i) = (k*dt/dxdx)*A(n,i-1) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,i) - [(udt/dx) - (k*dt/dxdx) ] * A(n,i+1)
    !
    ! k*dt/dxdx)*A(n,i-1) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,i) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,i+1) = A(n+1,i)
    !
    !(i= 1) =>  k*dt/dxdx)*A(n,0) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,1) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,2) = A(n+1,1)
    !(i= 2) =>  k*dt/dxdx)*A(n,1) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,2) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,3) = A(n+1,2)
    !(i= 3) =>  k*dt/dxdx)*A(n,2) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,3) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,4) = A(n+1,3)
    !(i= 4) =>  k*dt/dxdx)*A(n,3) + [1 + (udt/dx) - (2k*dt/dxdx)]*A(n,4) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,5) = A(n+1,4)
    !
    !(i= i) =>  k*dt/dxdx)*A(n,i-1) + [1 - (udt/dx) - (2k*dt/dxdx)]*A(n,i) + [(-udt/dx) + (k*dt/dxdx) ] * A(n,i+1) = A(n+1,i)
    !
    ! fronteria nao ciclica
    !   --                                                                             --         --    --      --      --
    !  | k*dt/dxdx)      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)  0  0   |    |A(n,0  ) |   | A(n+1,1) |
    !  | k*dt/dxdx)      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)         |    |A(n,1  ) |   | A(n+1,2) |
    !  | k*dt/dxdx)      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)         |  * |A(n,2  ) | = | A(n+1,3) |
    !  | k*dt/dxdx)      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)         |    |A(n,3  ) |   | A(n+1,4) |
    !
    !  | k*dt/dxdx)      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)         |    |A(n,i-1) |   | A(n+1,i) |
    !   --                                                                             --         --     --     --      --
    !
    ! fronteria ciclica
    !   --                                                                                                                                            --     --    --      --      --
    !  | [1 + (udt/dx) - (2k*dt/dxdx)]            [(-udt/dx) + (k*dt/dxdx)]                    0                       [k*dt/dxdx)]                     |    |A(n,0  ) |   | A(n+1,1) |
    !  | [ k*dt/dxdx)]                            [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)]          0                                |    |A(n,1  ) |   | A(n+1,2) |
    !  |      0                                   [k*dt/dxdx)]                      [1 + (udt/dx) - (2k*dt/dxdx)]     [(-udt/dx) + (k*dt/dxdx)]         |  * |A(n,2  ) | = | A(n+1,3) |
    !  |      0                                   0                                 [k*dt/dxdx) ]                     [1 + (udt/dx) - (2k*dt/dxdx)]     |    |A(n,3  ) |   | A(n+1,4) |
    !
    !  | [(-udt/dx) + (k*dt/dxdx)]                0                                 k*dt/dxdx)                        [1 + (udt/dx) - (2k*dt/dxdx)]     |    |A(n,i-1) |   | A(n+1,i) |
    !   --                                                                                                                                            --      --     --     --      --

    alpha = c*dt/dx
    beta  = kappa*dt/(dx*dx)

    do i = 2, nx-1
      ! indices com condicao periodica
      B(i,i+1) = -alpha + beta                             ! vizinho a esquerda
      B(i,i  ) = 1.0_real64 + alpha - 2.0_real64*beta     ! ponto central
      B(i,i-1) = beta                                     ! vizinho a direita
    end do
    ! indices com condicao periodica
    B(1,2 ) = beta                            ! vizinho a esquerda
    B(1,1 ) = 1.0_real64 + alpha - 2.0_real64*beta     ! ponto central
    B(1,nx) = -alpha + beta                                      ! vizinho a direita
      ! indices com condicao periodica
    B(nx, 1  ) = -alpha + beta                            ! vizinho a esquerda
    B(nx,nx  ) = 1.0_real64 + alpha - 2.0_real64*beta     ! ponto central
    B(nx,nx-1) = beta                                    ! vizinho a direita
    A=transpose(B) 
!PK   do i = 1, nx
!PK     ! indices com condicao periodica
!PK     j = i
!PK     A(i,j) = 1.0_real64 - alpha - 2.0_real64*beta     ! ponto central

!PK     j = i-1
!PK     if (j < 1) j = nx
!PK     A(i,j) =	alpha + beta				 ! vizinho a esquerda

!PK     j = i+1
!PK     if (j > nx) j = 1
!PK     A(i,j) =  beta					      ! vizinho a direita
!PK   end do

  end subroutine build_model_matrix
  !===========================

  !===========================
  subroutine compute_powers_of_A(A, Ap, nx, nmax)
    !! Ap(:,:,k) = A^k, para k=0..nmax
    real(real64), intent(in)  :: A(nx,nx)
    real(real64), intent(out) :: Ap(nx,nx,0:nmax)
    integer, intent(in)       :: nx, nmax

    integer :: k

    Ap(:,:,0) = identity_matrix(nx)
    Ap(:,:,1) = A
    do k = 2, nmax
      Ap(:,:,k) = matmul(Ap(:,:,k-1), A)
    end do

  end subroutine compute_powers_of_A
  !===========================

  !===========================
  subroutine integrate_truth(A, x0, traj, nx, nt)
    !! Integra x_{n+1} = A x_n
    real(real64), intent(in)    :: A(nx,nx)
    real(real64), intent(in)    :: x0(nx)
    real(real64), intent(out)   :: traj(nx,0:nt)
    integer, intent(in)         :: nx, nt

    integer :: k

    traj(:,0) = x0
    do k = 1, nt
      traj(:,k) = matmul(A, traj(:,k-1))
    end do

  end subroutine integrate_truth
  !===========================

  !===========================
  subroutine build_B_gaussian(B, Binv, x, nx, sigma_b, Lcorr)
    real(real64), intent(out) :: B(nx,nx), Binv(nx,nx)
    real(real64), intent(in)  :: x(nx)
    integer, intent(in)       :: nx
    real(real64), intent(in)  :: sigma_b, Lcorr

    integer :: i, j
    real(real64) :: d

    do i = 1, nx
      do j = 1, nx
        d = x(i) - x(j)
        B(i,j) = sigma_b**2 * exp( - (d*d) / (2.0_real64*Lcorr**2) )
      end do
    end do

    call invert_matrix(B, Binv, nx)

  end subroutine build_B_gaussian
  !===========================

  !===========================
  subroutine sample_background(x_true0, x0_bkg, B, nx)
    !! x0_bkg = x_true0 + b, com b ~ N(0,B)
    !!
    !! Aqui, para simplificar (e manter código compacto),
    !! vou gerar um ruído gaussiano i.i.d. com sigma_b
    !! em vez de amostrar exatamente de B.
    real(real64), intent(in)  :: x_true0(nx)
    real(real64), intent(out) :: x0_bkg(nx)
    real(real64), intent(in)  :: B(nx,nx)
    integer, intent(in)       :: nx

    integer :: i
    real(real64) :: rand1, rand2, z, sigma_b_local

    sigma_b_local = sqrt(B(1,1))   ! aproximação: var na diagonal ~ sigma_b^2

    do i = 1, nx
      call random_number(rand1)
      call random_number(rand2)
      if (rand1 <= 0.0_real64) rand1 = 1.0e-12_real64
      z = sqrt(-2.0_real64*log(rand1)) * cos(2.0_real64*pi*rand2)
      x0_bkg(i) = x_true0(i) + sigma_b_local * z
    end do

  end subroutine sample_background
  !===========================

  !===========================
  subroutine define_obs_points(obs_idx, nobs, nx)
    integer, intent(out) :: obs_idx(nobs)
    integer, intent(in)  :: nobs, nx

    integer :: i
    real(real64) :: step

    step = real(nx,real64) / real(nobs+1,real64)
    do i = 1, nobs
      obs_idx(i) = int(step*real(i,real64))
      if (obs_idx(i) < 1 ) obs_idx(i) = 1
      if (obs_idx(i) > nx) obs_idx(i) = nx
    end do

  end subroutine define_obs_points
  !===========================

  !===========================
  subroutine build_H_and_R(H, Ht, R, Rinv, obs_idx, nobs, nx, sigma_obs)
    real(real64), intent(out) :: H(nobs,nx), Ht(nx,nobs)
    real(real64), intent(out) :: R(nobs,nobs), Rinv(nobs,nobs)
    integer, intent(in)       :: obs_idx(nobs)
    integer, intent(in)       :: nobs, nx
    real(real64), intent(in)  :: sigma_obs

    integer :: i

    H = 0.0_real64
    do i = 1, nobs
      H(i, obs_idx(i)) = 1.0_real64
    end do
    Ht = transpose(H)

    R = 0.0_real64
    do i = 1, nobs
      R(i,i) = sigma_obs**2
    end do

    call invert_matrix(R, Rinv, nobs)

  end subroutine build_H_and_R
  !===========================

  !===========================
  subroutine generate_observations(truth_traj, y_obs, H, obs_idx, nobs, nx, sigma_obs)
    integer, intent(in)       :: nobs, nx
    real(real64), intent(in)  :: truth_traj(nx)
    real(real64), intent(in)  :: H(nobs,nx)
    integer, intent(in)       :: obs_idx(nobs)
    real(real64), intent(in)  :: sigma_obs
    real(real64), intent(out) :: y_obs(nobs)

    integer :: k, i
    real(real64) :: rand1, rand2, z

      do i = 1, nobs
        call random_number(rand1)
        call random_number(rand2)
        if (rand1 <= 0.0_real64) rand1 = 1.0e-12_real64
        z = sqrt(-2.0_real64*log(rand1)) * cos(2.0_real64*pi*rand2)
        y_obs(i) = truth_traj(obs_idx(i)) + sigma_obs*z
      end do

  end subroutine generate_observations
  !===========================

  !===========================
  subroutine compute_rmse(x1, x2, nx, rmse)
    real(real64), intent(in)  :: x1(nx), x2(nx)
    integer, intent(in)       :: nx
    real(real64), intent(out) :: rmse

    integer :: i
    real(real64) :: acc

    acc = 0.0_real64
    do i = 1, nx
      acc = acc + (x1(i)-x2(i))**2
    end do
    rmse = sqrt( acc / real(nx,real64) )

  end subroutine compute_rmse
  !===========================

  !===========================
  function identity_matrix(n) result(Id)
    integer, intent(in) :: n
    real(real64) :: Id(n,n)
    integer :: i

    Id = 0.0_real64
    do i = 1, n
      Id(i,i) = 1.0_real64
    end do
  end function identity_matrix
  !===========================

  !===========================
  subroutine invert_matrix(Ain, Ainv, n)
    !! Gauss-Jordan simples para n moderado
    real(real64), intent(in)  :: Ain(n,n)
    real(real64), intent(out) :: Ainv(n,n)
    integer, intent(in)       :: n

    real(real64), dimension(n,n) :: A
    real(real64), dimension(n,n) :: Id
    integer :: i, j, k, pivot
    real(real64) :: factor, maxval, tmp

    A = Ain
    Id = identity_matrix(n)

    do i = 1, n
      ! pivot parcial
      pivot = i
      maxval = abs(A(i,i))
      do k = i+1, n
        if (abs(A(k,i)) > maxval) then
          maxval = abs(A(k,i))
          pivot = k
        end if
      end do

      if (maxval < 1.0e-15_real64) then
        print *, 'ERRO: matriz quase singular na inversao, pivot=', i
        stop
      end if

      if (pivot /= i) then
        A([i,pivot],:)  = A([pivot,i],:)
        Id([i,pivot],:) = Id([pivot,i],:)
      end if

      tmp = A(i,i)
      A(i,:) = A(i,:) / tmp
      Id(i,:) = Id(i,:) / tmp

      do k = 1, n
        if (k /= i) then
          factor = A(k,i)
          A(k,:) = A(k,:) - factor*A(i,:)
          Id(k,:) = Id(k,:) - factor*Id(i,:)
        end if
      end do
    end do

    Ainv = Id

  end subroutine invert_matrix
  !===========================

  !===========================
  subroutine init_random_seed()
    integer :: i, nseed
    integer, allocatable :: seed(:)
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 1234567
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine init_random_seed
  !===========================

end program advec_diff_4dvar_equiv
