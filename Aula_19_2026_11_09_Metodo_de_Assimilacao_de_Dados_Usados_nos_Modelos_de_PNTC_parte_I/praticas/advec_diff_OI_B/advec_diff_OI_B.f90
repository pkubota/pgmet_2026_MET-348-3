!======================================================================
! Programa : advec_diff_OI_B.f90
! Descrição: Adveccao-difusao 1D + Assimilação por OI com B nao-diagonal
!            B: covariancia gaussiana espacial
! Padrão   : Fortran 2023
! Autor    : convertido por GPT-5 (exemplo)
!======================================================================

module advec_diff_mod
    use iso_fortran_env, only: real64
    implicit none
contains
    pure function step_advec_diff(c, nx, dx, dt, u, kappa) result(c_new)
        integer, intent(in) :: nx
        real(real64), intent(in) :: dx, dt, u, kappa
        real(real64), intent(in) :: c(nx)
        real(real64) :: c_new(nx)
        integer :: i
        real(real64) :: adv, diff

        c_new = c
        do concurrent (i = 2:nx-1)
            adv  = -u * (c(i) - c(i-1)) / dx
            diff =  kappa * (c(i+1) - 2.0_real64*c(i) + c(i-1)) / dx**2
            c_new(i) = c(i) + dt * (adv + diff)
        end do
    end function step_advec_diff
end module advec_diff_mod


module lin_alg_mod
    use iso_fortran_env, only: real64
    implicit none
contains
    function invert_matrix(A, n) result(Ainv)
        ! Gauss-Jordan inversion. Assumes A is non-singular.
        integer, intent(in) :: n
        real(real64), intent(in) :: A(n,n)
        real(real64) :: Ainv(n,n)
        real(real64), allocatable :: M(:,:)
        integer :: i, j, k
        real(real64) :: pivot, fac

        allocate(M(n,n))
        M = A
        Ainv = 0.0_real64
        do i = 1, n
            Ainv(i,i) = 1.0_real64
        end do

        do i = 1, n
            pivot = M(i,i)
            if (abs(pivot) < 1.0e-14_real64) then
                ! try to swap with a lower row
                do k = i+1, n
                    if (abs(M(k,i)) > abs(pivot)) then
                        M([i,k],:) = M([k,i],:)
                        Ainv([i,k],:) = Ainv([k,i],:)
                        pivot = M(i,i)
                        exit
                    end if
                end do
            end if
            if (abs(pivot) < 1.0e-14_real64) then
                stop "invert_matrix: singular matrix or pivot ~ 0"
            end if
            ! normalize row i
            M(i,:) = M(i,:) / pivot
            Ainv(i,:) = Ainv(i,:) / pivot

            ! eliminate other rows
            do j = 1, n
                if (j /= i) then
                    fac = M(j,i)
                    if (abs(fac) > 0.0_real64) then
                        M(j,:) = M(j,:) - fac * M(i,:)
                        Ainv(j,:) = Ainv(j,:) - fac * Ainv(i,:)
                    end if
                end if
            end do
        end do

        deallocate(M)
    end function invert_matrix


    pure function matvec_mul(A, x, n, m) result(y)
        ! y = A * x where A(n,m), x(m)
        integer, intent(in) :: n, m
        real(real64), intent(in) :: A(n,m), x(m)
        real(real64) :: y(n)
        integer :: i, j
        y = 0.0_real64
        do i = 1, n
            do j = 1, m
                y(i) = y(i) + A(i,j) * x(j)
            end do
        end do
    end function matvec_mul

end module lin_alg_mod


program advec_diff_OI_B
    USE Class_Fields, only : xa
    USE ModAdvection, only : AnaliticFunction
    USE Class_WritetoGrads, Only : SchemeWriteCtl,SchemeWriteData
    USE Model, only: Init,Run_Model,Finalize
    use iso_fortran_env, only: real64
    use advec_diff_mod
    use lin_alg_mod
    implicit none
    LOGICAL :: monotone=.true.

    ! Parameters
    integer, parameter      :: nx = 50, nt = 100
    real(real64), parameter :: dx = 1.0_real64 / nx, dt    = 0.001_real64
    real(real64), parameter :: u  = 1.0_real64     , kappa = 0.01_real64

    ! State arrays
    real(real64), allocatable :: x(:)
    real(real64), allocatable :: c_true(:), c_model(:), c_analysis(:)
!    real(real64), allocatable :: c0(:)

    ! Covariance / OI
    real(real64), allocatable :: B(:,:), B_sub(:,:), K(:,:), S(:,:), Sinv(:,:)

    real(real64), allocatable :: innov(:)
    real(real64), allocatable :: increment(:)
    real(real64), allocatable :: Tnew(:)
    ! Observations
    integer, PARAMETER   :: nobs=nx !observation number point 
    integer, allocatable :: obs_idx(:)
    real(real64), allocatable :: y_obs(:)

    real(real64) :: obs_error = 0.05_real64
    real(real64) :: model_error = 0.1_real64   ! sigma_b
    real(real64) :: Lcorr = 0.1_real64        ! comprimento de correlação
    real(real64) :: gain_diag

    ! temp & loop
    integer :: n,nn, i, j,n_step,irec,test
    real(real64) :: r1, r2,step_t,err
    REAL (KIND=real64) :: qmax
    REAL (KIND=real64) :: qmin

    ! ---------------------------
    ! Allocate
    ! ---------------------------
    allocate(x(nx),Tnew(nx), c_true( nx), c_model( nx), c_analysis(nx))
    !allocate(c0(nx))
    ! definir observacoes (mesmos índices do exemplo)
    allocate(obs_idx(nobs))
    allocate(y_obs(nobs))
    ! B e matrizes auxiliares
    allocate(B(nx, nx))
    allocate(B_sub(nx, nobs))
    allocate(K(nx, nobs))
    allocate(S(nobs, nobs))
    allocate(Sinv(nobs, nobs))
    allocate(innov(nobs))
    allocate(increment(nx))

    call Init(nx)
    ! ---------------------------
    ! Grade e condicao inicial
    ! ---------------------------
    do i = 1, nx
        x(i) = xa(i)   ! USE Class_Fields, only : xa real(i-1, real64) / real(nx-1, real64)
        !c0(i) = exp(-200.0_real64 * (x(i) - 0.3_real64)**2)
    end do
    !
    ! value of the observation
    !
    obs_idx=-999
    do i = 1, nx,1
         obs_idx(i) = i   ! Fortran 1-based
    end do
 

    irec=0
    step_t=0
    n_step=0 
    c_true (1:nx)  = AnaliticFunction(n_step)!c0(1:nx)                 !pluma de concentracao
    c_model(1:nx)  = Run_Model(step_t,n_step,real(nt,kind=real64),monotone)!0.8_real64 * c0(1:nx)    !pluma de concentracao

    test=SchemeWriteData(c_model (1:nx),irec)
    test=SchemeWriteData(c_true  (1:nx),irec)
    test=SchemeWriteData(c_model (1:nx),irec)

    ! ---------------------------
    ! Simulacao (verdade e modelo)
    ! ---------------------------
    err =0.0
    test=0
    do n = 2, nt
        c_true (1:nx) = AnaliticFunction(n_step)
        c_model(1:nx) = Run_Model(step_t,n_step,real(nt,kind=real64),monotone)

        !c_true (1:nx) = step_advec_diff(c_true (1:nx), nx, dx, dt, u, kappa)
        !c_model(1:nx) = step_advec_diff(c_model(1:nx), nx, dx, dt, u, kappa)
        !PK end do
        ! ---------------------------
        ! Construir B (covariancia gaussiana)
        ! B_ij = sigma_b^2 * exp( - (x_i - x_j)^2 / (2 L^2) )
        ! ---------------------------
        do i = 1, nx
           do j = 1, nx
              B(i,j) = model_error**2  *  exp( - (x(i)-x(j))**2 / (2.0_real64 * Lcorr**2) )
           end do
        end do

        ! ---------------------------
        ! Gerar observacoes (ruido Gaussiano) usando Box-Muller
        ! ---------------------------
        call random_seed()
        !PK do n = 1, nt
        do i = 1, nobs
            call random_number(r1)
            call random_number(r2)
            if (r1 <= 0.0_real64) r1 = 1.0e-12_real64
            ! Box-Muller ? z ~ N(0,1)
            !IF(obs_idx(i) >0)THEN
               !y_obs(i) = c_true(obs_idx(i)) + obs_error * sqrt(-2.0_real64*log(r1)) * cos(2.0_real64*acos(-1.0_real64)*r2)
               y_obs(i) = c_true(obs_idx(i)) + (c_model(obs_idx(i))-c_true (obs_idx(i)))/2.0

            !ELSE
            !   y_obs(i) = c_model(i)
            !end if
        end do
        !PK end do

        ! ---------------------------
        ! Preencher B_sub = B(:, obs_idx) (nx x nobs)
        ! ---------------------------
        do i = 1, nx
           do j = 1, nobs
              IF(obs_idx(j) >0)THEN
                 B_sub(i,j) = B(i, obs_idx(j))
              ELSE 
                B_sub(i,j) = 0.0
              END IF
           end do
        end do

        ! ---------------------------
        ! Calcular S = H B H^T + R  (nobs x nobs)
        ! Como H seleciona indices, HBHT = B(obs_idx, obs_idx)
        ! R = obs_error^2 * I
        ! ---------------------------
        do i = 1, nobs
           do j = 1, nobs
              IF(obs_idx(i) >0 .and.obs_idx(j) >0 )THEN
                 S(i,j) = B(obs_idx(i), obs_idx(j))
              ELSE 
                 S(i,j) = 0.0
              END IF

           end do
        end do

        do i = 1, nobs
           S(i,i) = S(i,i) + obs_error**2
        end do

        ! Inverter S 
        Sinv(1:nobs, 1:nobs) = invert_matrix(S(1:nobs, 1:nobs), nobs)

        ! ---------------------------
        ! Calcular ganho K = B_sub * Sinv   (nx x nobs) * (nobs x nobs) -> (nx x nobs)
        ! ---------------------------
        do i = 1, nx
           do j = 1, nobs
              K(i,j) = 0.0_real64
              do nn = 1, nobs
                 K(i,j) = K(i,j) + B_sub(i,nn) * Sinv(nn,j)
              end do
           end do
        end do

        ! ---------------------------
        ! Aplicar OI em cada passo temporal:
        ! xa = xb + K (y - H xb)
        ! H xb e apenas xb nos indices observados
        ! ---------------------------
        c_analysis(1:nx) = c_model(1:nx)
        !PK do n = 1, nt

        do i = 1, nobs
              IF(obs_idx(i) >0)THEN
                 innov(i) = y_obs(i) - c_model(obs_idx(i))
              ELSE 
                 innov(i) = 0.0
              END IF
        end do

        ! increment = K * innov -> vetor nx

        increment(1:nx) = matvec_mul(K(1:nx, 1:nobs), innov(1:nobs), nx, nobs)
        ! atualizar análise
         do i = 1, nx
            c_analysis(i) = c_model(i) + increment(i)
         end do
        if( monotone )THEN
         Tnew(1)=c_analysis(1)
         qmax = max(c_analysis(nx),c_analysis(1))
         qmin = min(c_analysis(nx),c_analysis(1))
         Tnew(1) = max(qmin,min(qmax,c_analysis(1)))
         do i = 2, nx
            qmax = max(c_analysis(i-1),c_analysis(i))
            qmin = min(c_analysis(i-1),c_analysis(i))
            Tnew(i) = max(qmin,min(qmax,c_analysis(i)))
         end do
         c_analysis=Tnew
        else
           do i = 1, nx
               c_analysis(i) = c_model(i) + increment(i)
           end do
        end if
        test=SchemeWriteData(c_model    (1:nx),irec)
        test=SchemeWriteData(c_true     (1:nx),irec)
        test=SchemeWriteData(c_analysis (1:nx),irec)

    end do
    test=SchemeWriteCtl(n_step)


    ! ---------------------------
    ! Saida: imprimir alguns valores e salvar arquivo para plot
    ! ---------------------------
    write(*,'(A)') '=========================================================='
    write(*,'(A)') '  OI com B nao-diagonal (covariancia gaussiana espacial)'
    write(*,'(A)') '=========================================================='
    write(*,'(A,F8.4)') 'sigma_b (model_error) = ', model_error
    write(*,'(A,F8.4)') 'Lcorr = ', Lcorr
    write(*,'(A,F8.4)') 'obs_error = ', obs_error
    write(*,'(A)') 'Inovacao e analise final nos pontos observados:'
    do i = 1, nobs
        IF(obs_idx(i) >0)THEN
        write(*,'(A,I3, A, F8.4, A, F8.4)') ' idx=', obs_idx(i), &
            '   obs=', y_obs(i), '   analise=', c_analysis(obs_idx(i))
        end if
    end do
    write(*,'(A)') 'Arquivo "oi_B_output.dat" gerado para plot externo.'
    write(*,'(A)') '=========================================================='

    open(unit=11, file='oi_B_output.dat', status='replace')
    write(11,'(A)') '# x  c_true  c_model  c_analysis'
    do i = 1, nx
        write(11,'(4F12.6)') x(i), c_true(i), c_model(i), c_analysis(i)
    end do
    close(11)

    ! liberar
    deallocate(B, B_sub, K, S, Sinv)
    deallocate(x, c_true, c_model, c_analysis, y_obs, obs_idx)
    !deallocate(c0)
    deallocate(innov)
    deallocate(increment)

end program advec_diff_OI_B
