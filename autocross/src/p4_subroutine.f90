! Subroutine interface to PearsonT4. See comments in pearsont3.f90
!

subroutine p4_subroutine(n, t, x, y, confidence_level, b1_in, b2_in, &
  n_lambda_in, corr, ci, taux, tauy)
  use result1
  use pearsont3_module
  use data1
  use data2
  use random
  use time
  use setting, only: n1
  use parameters, only: alpha
  implicit none
  integer, intent(in) :: n,b1_in,b2_in,n_lambda_in
  real(dp), intent(in) :: t(n), x(n), y(n), confidence_level
  real(dp), intent(out) :: corr, ci(2), taux, tauy
  logical :: allocate0_here
  logical :: profile
  character(len=32) :: profile_env
  character(len=32) :: seed_env
  character(len=32) :: b1_env,b2_env,b3_env
  integer :: base_seed, seed_status
  integer :: b1_read,b2_read,b3_read,read_status
  integer, save :: p4_call_count = 0
  real :: t_phase, t_total
  call get_environment_variable("P4_PROFILE", profile_env)
  profile = len_trim(profile_env) > 0
  if (profile) call cpu_time(t_total)
  write(*,*) 'p4_subroutine: start, n =', n
  call flush(6)
  ! The public interface uses the ordinary two-sided significance level
  ! (0.05 for a 95% CI). PearsonT4's internal alpha is one tail.
  alpha = 0.5_dp * confidence_level
  ! PearsonT4 does not use the observed-data T3 B1/B2 bootstrap.
  call set_bootstrap_parameters(1,1,n_lambda_in)
  call set_pearson_t4_null_calibration(.true.)
  b1_read = 1000
  b2_read = 100
  b3_read = 100
  call get_environment_variable("P4_B1", b1_env)
  if (len_trim(b1_env) > 0) then
    read(b1_env, *, iostat=read_status) b1_read
    if (read_status /= 0) b1_read = 1000
  end if
  call get_environment_variable("P4_B2", b2_env)
  if (len_trim(b2_env) > 0) then
    read(b2_env, *, iostat=read_status) b2_read
    if (read_status /= 0) b2_read = 100
  end if
  call get_environment_variable("P4_B3", b3_env)
  if (len_trim(b3_env) > 0) then
    read(b3_env, *, iostat=read_status) b3_read
    if (read_status /= 0) b3_read = 100
  end if
  call set_pearson_t4_calibration_parameters(b1_read,b2_read,b3_read)
  call get_environment_variable("P4_SEED", seed_env)
  if (len_trim(seed_env) == 0) call get_environment_variable("P3_SEED", seed_env)
  if (len_trim(seed_env) > 0) then
    read(seed_env, *, iostat=seed_status) base_seed
    if (seed_status == 0) then
      call ranseed(base_seed + p4_call_count)
      call seed_randompack(base_seed + p4_call_count)
      p4_call_count = p4_call_count + 1
    else
      call ranseed()
      call seed_randompack(1234567)
    end if
  else
    call ranseed()
    call seed_randompack(1234567)
  end if
  n1 = n
  allocate0_here = .not.allocated(t1)
  write(*,*) 'p4_subroutine: allocate0_here =', allocate0_here
  call flush(6)
  if (allocate0_here) then
    write(*,*) 'p4_subroutine: calling allocate0'
    call flush(6)
    if (profile) call cpu_time(t_phase)
    call allocate0  ! t1, x1, y1
    if (profile) call print_profile("allocate0", t_phase)
  end if
  t1 = t
  x1 = x
  y1 = y
  n1 = n
  write(*,*) 'p4_subroutine: calling init1a'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call init1a           ! n2=n1
  if (profile) call print_profile("init1a", t_phase)
  write(*,*) 'p4_subroutine: calling allocate1'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call allocate1              ! t2, x2, y2, x3, y3, x3_resample1, y3_resample1,
  !                             x3_resample2, Dy3_resample2
  if (profile) call print_profile("allocate1", t_phase)
  write(*,*) 'p4_subroutine: calling allocate_resample_data'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call allocate_resample_data
  if (profile) call print_profile("allocate_resample_data", t_phase)
  write(*,*) 'p4_subroutine: calling calc_t_inv_lambda'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call calc_t_inv_lambda      ! calculates percentage point tv(lambda) over a
  !                             lambda grid (Calibrated CI)
  if (profile) call print_profile("calc_t_inv_lambda", t_phase)
  write(*,*) 'p4_subroutine: calling init1b'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call init1b                 ! t2, x2, y2, x3, y3, x3_resample1,
  !                             y3_resample1,
  !                             x3_resample2, y3_resample2
  if (profile) call print_profile("init1b", t_phase)
  write(*,*) 'p4_subroutine: calling r_est'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call r_est       ! detrends (x2->x3, y2->y3), estimates r(x3, y3)
  if (profile) call print_profile("r_est", t_phase)
  write(*,*) 'p4_subroutine: calling tauest'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call tauest      ! estimates persistence times taux3, tauy3 and rhox3 and
  !                  rhoy3)
  if (profile) call print_profile("tauest", t_phase)
  write(*,*) 'p4_subroutine: calling chsett4'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call chsett4     ! changes setting: l_mbb , block length
  if (profile) call print_profile("chsett4", t_phase)
  write(*,*) 'p4_subroutine: calling confidence'
  call flush(6)
  if (profile) call cpu_time(t_phase)
  call confidence  ! estimates [r_low; r_upp]
  if (profile) call print_profile("confidence", t_phase)
  call set_pearson_t4_null_calibration(.false.)
  corr = r
  ci = [r_low, r_upp]
  taux = taux3
  tauy = tauy3
  write(*,*) 'p4_subroutine: before deallocation, allocate0_here =', allocate0_here
  call flush(6)
  if (allocate0_here) then
    write(*,*) 'p4_subroutine: deallocating temporary arrays'
    call flush(6)
    if (profile) call cpu_time(t_phase)
    call deallocate0
    call deallocate_resample_data
    call deallocate1
    if (profile) call print_profile("deallocate", t_phase)
  end if
  write(*,*) 'p4_subroutine: end'
  call flush(6)
  if (profile) call print_profile("total", t_total)
contains
  subroutine print_profile(name, started)
    character(len=*), intent(in) :: name
    real, intent(in) :: started
    real :: finished
    call cpu_time(finished)
    write(*,'("P4_PROFILE ",A,": ",F0.6," s")') trim(name), finished - started
  end subroutine print_profile
end subroutine p4_subroutine
