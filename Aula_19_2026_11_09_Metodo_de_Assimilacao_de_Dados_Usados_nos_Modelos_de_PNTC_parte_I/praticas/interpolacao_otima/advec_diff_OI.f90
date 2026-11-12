!======================================================================
! Programa : advec_diff_OI.f90
! Descrição: Modelo 1D de adveccao-difusao com assimilacao (OI)
! Autor    : Paulo Kubota 
! Padrão   : Fortran 2023
!======================================================================

module advec_diff_mod
    use iso_fortran_env, only: real64
    implicit none
contains
    pure function step_advec_diff(c, nx, dx, dt, u, kappa) result(c_new)
        ! Passo de tempo explicito da equacao de adveccao-difusao 1D
        integer, intent(in) :: nx
        real(real64), intent(in) :: dx, dt, u, kappa
        real(real64), intent(in) :: c(nx)
        real(real64) :: c_new(nx)
        integer :: i
        real(real64) :: adv, diff

        c_new = c
        do concurrent (i = 2:nx-1)
            adv  =  -u * (c(i) - c(i-1)) / dx
            diff =  kappa * (c(i+1) - 2.0_real64*c(i) + c(i-1)) / dx**2
            c_new(i) = c(i) + dt * (adv + diff)
        end do
    end function step_advec_diff
end module advec_diff_mod


program advec_diff_OI
    use iso_fortran_env, only: real64
    use advec_diff_mod
    implicit none

    ! ---------------------------
    ! Parametros do modelo
    ! ---------------------------
    integer, parameter :: nx = 50, nt = 50
    real(real64), parameter :: dx = 1.0_real64 / nx, dt = 0.001_real64
    real(real64), parameter :: u = 1.0_real64, kappa = 0.01_real64

    real(real64), allocatable :: x(:)
    real(real64), allocatable :: c_true(:,:), c_model(:,:), c_analysis(:,:)
    real(real64), allocatable :: c0(:)

    integer :: n, i
    integer, dimension(3) :: obs_idx = [11, 26, 41]  ! Fortran e 1-based
    real(real64) :: obs_error = 0.05_real64
    real(real64) :: model_error = 0.1_real64
    real(real64) :: gain
    real(real64), allocatable :: y_obs(:,:)
    real(real64) :: r

    ! ---------------------------
    ! Inicializacao
    ! ---------------------------
    allocate(x(nx), c_true(nt, nx), c_model(nt, nx), c_analysis(nt, nx))
    allocate(c0(nx), y_obs(nt, size(obs_idx)))

    ! Grade e condicao inicial (pulso gaussiano)
    do i = 1, nx
        x(i) = real(i-1, real64) / real(nx-1, real64)
        c0(i) = exp(-200.0_real64 * (x(i) - 0.3_real64)**2)
    end do

    c_true(1,:) = c0
    c_model(1,:) = 0.8_real64 * c0

    ! ---------------------------
    ! Simulacao da "verdade"
    ! ---------------------------
    do n = 2, nt
        c_true(n,:) = step_advec_diff(c_true(n-1,:), nx, dx, dt, u, kappa)
    end do

    ! ---------------------------
    ! Simulacao do "modelo"
    ! ---------------------------
    do n = 2, nt
        c_model(n,:) = step_advec_diff(c_model(n-1,:), nx, dx, dt, u, kappa)
    end do

    ! ---------------------------
    ! Geracao de observacoes (com ruido)
    ! ---------------------------
    call random_seed()  ! inicializa gerador
    do n = 1, nt
        do i = 1, size(obs_idx)
            call random_number(r)
            y_obs(n, i) = c_true(n, obs_idx(i)) + obs_error * (2.0_real64*r - 1.0_real64)
        end do
    end do

    ! ---------------------------
    ! Interpolacao Otima (OI)
    ! ---------------------------
    gain = model_error / (model_error + obs_error)
    c_analysis = c_model

    do n = 1, nt
        do i = 1, size(obs_idx)
            c_analysis(n, obs_idx(i)) = (1.0_real64 - gain)*c_model(n, obs_idx(i)) + gain*y_obs(n, i)
        end do
    end do

    ! ---------------------------
    ! Saída final (ultimo passo)
    ! ---------------------------
    write(*,'(A)') '==============================================='
    write(*,'(A)') ' Interpolacao Otima (OI) - Adveccao-Difusao 1D'
    write(*,'(A)') '==============================================='
    write(*,'(A,F6.3)') 'Ganho = ', gain
    write(*,'(A)') 'Pontos observados e analise final:'
    do i = 1, size(obs_idx)
        write(*,'(A,I3,2(A,F6.3))') ' idx=', obs_idx(i), '  obs=', y_obs(nt,i), &
                                     '  analise=', c_analysis(nt, obs_idx(i))
    end do
    write(*,'(A)') '==============================================='

    ! Exportar resultados (para plot posterior em Python ou Gnuplot)
    open(unit=10, file='oi_output.dat', status='replace')
    write(10,'(A)') '# x  c_true  c_model  c_analysis'
    do i = 1, nx
        write(10,'(4F12.6)') x(i), c_true(nt,i), c_model(nt,i), c_analysis(nt,i)
    end do
    close(10)
    write(*,*) 'Arquivo "oi_output.dat" gerado para plot externo.'
end program advec_diff_OI
