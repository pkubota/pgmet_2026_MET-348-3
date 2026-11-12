!===========================================================
! Programa: interpolacao_otima.f90
! Autor: Paulo Kubota 
!===========================================================

module class_interpolacao_otima_mod
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: optimal_interpolation
contains

    pure function optimal_interpolation(x_model, y_model, x_obs, y_obs, model_error, obs_error) result(y_analysis)
        ! Interpolacao otima simples (1D)
        real(real64), intent(in) :: x_model(:), y_model(:)
        real(real64), intent(in) :: x_obs(:), y_obs(:)
        real(real64), intent(in) :: model_error, obs_error
        real(real64) :: y_analysis(size(x_obs))

        real(real64) :: gain

        ! Ganha de interpolacao (simplificada)
        gain = model_error / (model_error + obs_error)

        ! Corrige o valor no ponto x = 2 ? índice 3 (pois Fortran e 1-based)
        y_analysis = (1.0_real64 - gain) * y_model(3) + gain * y_obs

    end function optimal_interpolation

end module class_interpolacao_otima_mod


program interpolacao_otima
    use iso_fortran_env, only: real64
    use class_interpolacao_otima_mod
    implicit none

    real(real64), dimension(6) :: x_model, y_model
    real(real64), dimension(1) :: x_obs, y_obs, y_analysis
    real(real64) :: model_error, obs_error

    ! Dados do modelo
    x_model = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    y_model = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0]

    ! Observacao
    x_obs = [2.0]
    y_obs = [15.5]

    ! Erros
    model_error = 1.0
    obs_error   = 0.5

    ! Aplicar Interpolacao Otima
    y_analysis = optimal_interpolation(x_model, y_model, x_obs, y_obs, model_error, obs_error)

    ! Saída dos resultados
    write(*,'(A)') '=========================================='
    write(*,'(A)') '      INTERPOLACAO OTIMA (OI) - FORTRAN   '
    write(*,'(A)') '=========================================='
    write(*,'(A,F6.2)') 'Ponto observado (x = ', x_obs(1)
    write(*,'(A,F6.2)') 'Temperatura observada: ', y_obs(1)
    write(*,'(A,F6.2)') 'Temperatura do modelo: ', y_model(3)
    write(*,'(A,F6.2)') 'Temperatura apos assimilacao: ', y_analysis(1)
    write(*,'(A)') '=========================================='

end program interpolacao_otima
