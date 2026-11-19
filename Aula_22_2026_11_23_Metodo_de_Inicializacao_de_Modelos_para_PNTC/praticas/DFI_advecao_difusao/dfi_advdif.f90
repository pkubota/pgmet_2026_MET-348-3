!===============================================================================
! Modelo 1D de adveccao-difusao com Inicializacao por Filtro Digital (DFI)
! no estilo do HIRLAM (Lynch e Huang, 1992), com filtro de Dolph-Chebyshev.
!
!     dc/dt + U dc/dx = nu d2c/dx2        (dominio periodico)
!
! Esquema: leapfrog centrado no espaco; difusao defasada (lagged); primeiro
! passo por Runge-Kutta de ponto medio. O leapfrog e reversivel no tempo
! (dt < 0), o que permite a integracao para tras do estagio 1 da DFI.
!
! DFI em dois estagios (como no HIRLAM):
!   estagio 1: de t=0 a t=-TS, para tras, SEM difusao (processo irreversivel),
!              acumulando sum_k h_k c(t_k) -> cbar valido em t=-TS/2
!   estagio 2: de t=-TS/2 a t=+TS/2, para frente, a partir de cbar,
!              acumulando os mesmos pesos -> campo inicializado valido em t=0
! A variante "diab" liga a difusao no estagio 2.
!
! Compilar:  gfortran -O2 -o dfi_advdif dfi_advdif.f90
! Executar:  ./dfi_advdif [diretorio_de_saida]
!===============================================================================
program dfi_advdif
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  ! ------------------------------------------------------------- parametros
  integer,  parameter :: nx  = 400
  real(dp), parameter :: xlen = 400.0e3_dp        ! comprimento do dominio (m)
  real(dp), parameter :: u   = 10.0_dp           ! velocidade (m/s)
  real(dp), parameter :: nu  = 100.0_dp          ! difusividade (m2/s)
  real(dp), parameter :: dt  = 50.0_dp           ! passo de tempo (s); Courant=0.5
  real(dp), parameter :: tau_s = 900.0_dp        ! periodo de corte (s)
  integer,  parameter :: mm  = 10                ! N = 2M+1 pesos; TS = 2 M dt
  real(dp), parameter :: x0c = 120.0e3_dp        ! centro da gaussiana (m)
  real(dp), parameter :: sig0 = 20.0e3_dp        ! largura da gaussiana (m)
  integer,  parameter :: nfc = 216               ! passos da previsao (3 h)
  real(dp), parameter :: noise_amp = 0.008_dp    ! amplitude por modo do ruido

  real(dp) :: dx, theta_s, x0, r
  real(dp) :: x(nx), ctrue(nx), cnil(nx), cdfi(nx), cdfd(nx)
  real(dp) :: f_nil(nx), f_dfi(nx), f_dfd(nx), ctend(nx)
  real(dp) :: h(0:2*mm)
  real(dp) :: r_nil(0:nfc), r_dfi(0:nfc), r_dfd(0:nfc)
  character(len=256) :: outdir
  integer :: i, k, narg

  narg = command_argument_count()
  outdir = 'saida_f90'
  if (narg >= 1) call get_command_argument(1, outdir)
  call execute_command_line('mkdir -p '//trim(outdir))

  dx = xlen / real(nx, dp)

  do i = 1, nx
     x(i) = (real(i - 1, dp) + 0.5_dp) * dx
  end do

  theta_s = 2.0_dp * pi * dt / tau_s

  call dolph(mm, theta_s, h, x0, r)

  write (*, '(A,F6.3,A,F8.5,A,I3,A,I3,A,F6.0,A)') 'Courant=', u*dt/dx, '  theta_s=', theta_s, &
       '  M=', mm, '  N=', 2*mm+1, '  TS=', 2.0_dp*mm*dt, ' s'
  write (*, '(A,F9.6,A,F8.5,A,F6.2,A,F15.12,A,ES11.4)') 'x0=', x0, '  r=', r, ' (', &
       20.0_dp*log10(r), ' dB)  soma(h)=', sum(h), '  h_min=', minval(h)

  call analitica(x, 0.0_dp, ctrue)
  call ruido(x, ctend)

  cnil = ctrue + ctend                           ! analise "nao inicializada"

  call dfi(cnil, h, .false., cdfi)               ! DFI adiabatica (HIRLAM)
  call dfi(cnil, h, .true.,  cdfd)               ! variante com difusao no estagio 2

  call previsao(cnil, r_nil, f_nil)
  call previsao(cdfi, r_dfi, f_dfi)
  call previsao(cdfd, r_dfd, f_dfd)

  write (*, '(A,3F10.5)') 'RMSE em t=0 (NIL, DFI, DFI-diab):', rms(cnil-ctrue), rms(cdfi-ctrue), rms(cdfd-ctrue)

  do k = 0, nfc, 36
     if (k == 36 .or. k == 108 .or. k == 0 .or. k == 216) then
        write (*, '(A,F4.1,A,3F10.5)') 'RMSE t=', real(k,dp)*dt/3600.0_dp, ' h:', r_nil(k), r_dfi(k), r_dfd(k)
     end if
  end do

  ! ------------------------------------------------------------- saidas
  open (10, file=trim(outdir)//'/campos_t0.dat', status='replace')
  write (10, '(A)') '# x ctrue cnil cdfi cdfi_diab'
  do i = 1, nx
     write (10, '(5ES26.15E3)') x(i), ctrue(i), cnil(i), cdfi(i), cdfd(i)
  end do
  close (10)


  open (10, file=trim(outdir)//'/rmse.dat', status='replace')
  write (10, '(A)') '# t_s rmse_nil rmse_dfi rmse_dfi_diab'
  do k = 0, nfc
     write (10, '(4ES26.15E3)') real(k,dp)*dt, r_nil(k), r_dfi(k), r_dfd(k)
  end do
  close (10)

  open (10, file=trim(outdir)//'/pesos.dat', status='replace')
  write (10, '(A)') '# k h_k'
  do k = -mm, mm
     write (10, '(I4,ES26.15E3)') k, h(k+mm)
  end do
  close (10)

  call analitica(x, real(nfc,dp)*dt, ctend)
  open (10, file=trim(outdir)//'/campos_prev.dat', status='replace')
  write (10, '(A)') '# x c_nil c_dfi c_dfi_diab c_true (fim da previsao)'
  do i = 1, nx
     write (10, '(5ES26.15E3)') x(i), f_nil(i), f_dfi(i), f_dfd(i), ctend(i)
  end do
  close (10)

contains

  !---------------------------------------------------------- utilitarios
  pure function rms(a) result(v)
    real(dp), intent(in) :: a(:)
    real(dp) :: v
    v = sqrt(sum(a**2) / real(size(a), dp))
  end function rms

  !---------------------------------------------- polinomio de Chebyshev T_n(x)
  pure function cheb(n, xx) result(t)
    integer,  intent(in) :: n
    real(dp), intent(in) :: xx
    real(dp) :: t
    if (abs(xx) <= 1.0_dp) then
       t = cos(real(n, dp) * acos(xx))
    else
       t = cosh(real(n, dp) * acosh(xx))
    end if
  end function cheb

  !---------------------------------------- pesos do filtro de Dolph-Chebyshev
  ! h(s), s = 0..2M corresponde a k = s - M (pesos simetricos)
  subroutine dolph(m, ths, hh, xx0, rr)
    integer,  intent(in)  :: m
    real(dp), intent(in)  :: ths
    real(dp), intent(out) :: hh(0:2*m), xx0, rr
    integer  :: nn, kk, mi
    real(dp) :: thm, sm
    nn  = 2 * m + 1
    xx0 = 1.0_dp / cos(ths / 2.0_dp)
    rr  = 1.0_dp / cosh(2.0_dp * real(m, dp) * acosh(xx0))
    do kk = -m, m
       sm = 0.0_dp
       do mi = 1, m
          thm = 2.0_dp * pi * real(mi, dp) / real(nn, dp)
          sm  = sm + cheb(2*m, xx0 * cos(thm / 2.0_dp)) * cos(real(kk, dp) * thm)
       end do
       hh(kk + m) = (1.0_dp + 2.0_dp * rr * sm) / real(nn, dp)
    end do
  end subroutine dolph

  !-------------------------------- solucao analitica (imagens periodicas)
  subroutine analitica(xg, t, c)
    real(dp), intent(in)  :: xg(:), t
    real(dp), intent(out) :: c(:)
    real(dp) :: s2, amp, xc
    integer  :: j
    s2  = sig0**2 + 2.0_dp * nu * t
    amp = sig0 / sqrt(s2)
    c = 0.0_dp
    do j = -1, 1
       xc = x0c + u * t + real(j, dp) * xlen
       c = c + amp * exp(-0.5_dp * (xg - xc)**2 / s2)
    end do
  end subroutine analitica

  !------------------------------ ruido deterministico (ondas de 8dx a 3dx)
  subroutine ruido(xg, n)
    real(dp), intent(in)  :: xg(:)
    real(dp), intent(out) :: n(:)
    integer  :: m
    real(dp) :: ph
    n = 0.0_dp
    do m = nx/8, nx/3  ! 50..133 para nx=400
       ph = 2.0_dp * pi * modulo(real(m, dp) * 0.6180339887498949_dp, 1.0_dp)
       n = n + noise_amp * cos(2.0_dp * pi * real(m, dp) * xg / xlen + ph)
    end do
  end subroutine ruido

  !---------------------------------------------------- tendencia espacial
  subroutine tend(c, adv, dif, f)
    real(dp), intent(in)  :: c(nx)
    logical,  intent(in)  :: adv, dif
    real(dp), intent(out) :: f(nx)
    integer :: j, jp, jm
    f = 0.0_dp
    do j = 1, nx
       jp = j + 1; if (jp > nx) jp = 1
       jm = j - 1; if (jm < 1)  jm = nx
       if (adv) f(j) = f(j) - u * (c(jp) - c(jm)) / (2.0_dp * dx)
       if (dif) f(j) = f(j) + nu * (c(jp) - 2.0_dp * c(j) + c(jm)) / dx**2
    end do
  end subroutine tend

  !------------------------------------------------------ primeiro passo (RK2)
  subroutine primeiro_passo(c, dtl, dif, c1)
    real(dp), intent(in)  :: c(nx), dtl
    logical,  intent(in)  :: dif
    real(dp), intent(out) :: c1(nx)
    real(dp) :: f(nx), cs(nx)
    call tend(c, .true., dif, f)
    cs = c + 0.5_dp * dtl * f
    call tend(cs, .true., dif, f)
    c1 = c + dtl * f
  end subroutine primeiro_passo

  !--------------------------------------------------- integracao leapfrog
  ! dtl pode ser negativo. Se h e acc estao presentes, acumula h(s)*c(s)
  ! (s = 0..ns). Se rmse esta presente, guarda RMSE(s) contra a analitica.
  subroutine integra(c0, dtl, ns, dif, cfinal, h, acc, rmse)
    real(dp), intent(in)    :: c0(nx), dtl
    integer,  intent(in)    :: ns
    logical,  intent(in)    :: dif
    real(dp), intent(out)   :: cfinal(nx)
    real(dp), intent(in),    optional :: h(0:)
    real(dp), intent(inout), optional :: acc(nx)
    real(dp), intent(inout), optional :: rmse(0:)
    real(dp) :: cold(nx), cnow(nx), cnew(nx), fa(nx), fd(nx)
    integer  :: s

    cold = c0
    call passo_diag(0, cold, dtl, h, acc, rmse)
    cnow = cold
    if (ns >= 1) then
       call primeiro_passo(cold, dtl, dif, cnow)
       call passo_diag(1, cnow, dtl, h, acc, rmse)
    end if
    do s = 2, ns
       call tend(cnow, .true., .false., fa)
       cnew = cold + 2.0_dp * dtl * fa
       if (dif) then
          call tend(cold, .false., .true., fd)
          cnew = cnew + 2.0_dp * dtl * fd
       end if
       cold = cnow
       cnow = cnew
       call passo_diag(s, cnow, dtl, h, acc, rmse)
    end do
    cfinal = cnow

  end subroutine integra

  !------------------- acumulo de pesos e/ou RMSE em cada passo da integracao
  subroutine passo_diag(ss, c, dtl, h, acc, rmse)
    integer,  intent(in) :: ss
    real(dp), intent(in) :: c(nx), dtl
    real(dp), intent(in),    optional :: h(0:)
    real(dp), intent(inout), optional :: acc(nx)
    real(dp), intent(inout), optional :: rmse(0:)
    real(dp) :: ct(nx)
    if (present(h) .and. present(acc)) acc = acc + h(ss) * c
    if (present(rmse)) then
       call analitica(x, real(ss, dp) * dtl, ct)
       rmse(ss) = rms(c - ct)
    end if
  end subroutine passo_diag

  !----------------------------------------------------------------- DFI
  subroutine dfi(c0, hh, diab, cinit)
    real(dp), intent(in)  :: c0(nx), hh(0:)
    logical,  intent(in)  :: diab
    real(dp), intent(out) :: cinit(nx)
    real(dp) :: acc(nx), cbar(nx), cdum(nx)
    integer  :: ns
    ns = 2 * mm
    ! estagio 1: para tras (dt<0), sem difusao
    acc = 0.0_dp
    call integra(c0, -dt, ns, .false., cdum, h=hh, acc=acc)
    cbar = acc                          ! valido em t = -TS/2
    ! estagio 2: para frente, a partir de cbar
    acc = 0.0_dp
    call integra(cbar, dt, ns, diab, cdum, h=hh, acc=acc)
    cinit = acc                         ! valido em t = 0
  end subroutine dfi

  !------------------------------------------------------------- previsao
  subroutine previsao(c0, rmse, cfim)
    real(dp), intent(in)  :: c0(nx)
    real(dp), intent(out) :: rmse(0:nfc), cfim(nx)
    call integra(c0, dt, nfc, .true., cfim, rmse=rmse)
  end subroutine previsao

end program dfi_advdif
