! Internal C ABI for the Python extension. The legacy R entry point is unchanged.
subroutine p3_python(n, t, x, y, nlevel, levels, nb1, nb2, nl, block, linear, draw, context, &
    corr, ci, tau, actual_block, se_z, selected, cover, z1, se2, grid, persistence, degenerate, zero_inner, status) bind(C)
  use iso_c_binding, only: c_int, c_int64_t, c_double, c_funptr, c_f_procpointer
  use pearsont3_module
  use parameters, only: alpha
  use data2
  use inv_tlambda, only: t_inv_lambda
  use resample_data
  use result1
  use setting, only: n2, dtrtype, l_mbb
  use minls_module, only: tauest_x
  implicit none
  integer(c_int), value :: n, nlevel, nb1, nb2, nl, block, linear
  type(c_funptr), value :: draw
  type(c_ptr), value :: context
  real(c_double), intent(in) :: t(n), x(n), y(n), levels(nlevel)
  real(c_double), intent(out) :: corr, ci(2,nlevel), tau(2), se_z, selected(nlevel)
  real(c_double), intent(out) :: cover(nl), z1(nb1), se2(nb1), grid(nl)
  integer(c_int), intent(out) :: actual_block, persistence(2), zero_inner, status
  integer(c_int64_t), intent(out) :: degenerate
  integer :: ierr, i, k
  real(dp) :: one_tail, target, z, half_width

  api_call = .true.
  api_status = 0
  api_degenerate = 0
  api_rng_context = context
  call c_f_procpointer(draw, api_draw)
  call set_bootstrap_parameters(nb1, nb2, nl)
  call set_pearson_t4_null_calibration(.false.)
  profile_detail_enabled = .false.
  alpha = (1.0_dp - levels(1)) / 2.0_dp
  n2 = n
  dtrtype = 'm'
  if (linear /= 0) dtrtype = 'l'

  status = 1 ! allocation failure
  call allocate1(ierr)
  if (ierr /= 0) goto 900
  call allocate_resample_data(ierr)
  if (ierr /= 0) goto 900
  allocate(rp_starts(n), stat=ierr)
  if (ierr /= 0) goto 900
  t2(:) = t
  x2(:) = x
  y2(:) = y
  call r_est
  status = 3 ! numerical failure
  if (.not. all(ieee_is_finite(x3)) .or. .not. all(ieee_is_finite(y3))) goto 900
  status = 5 ! constant after detrending
  if (all(x3 == x3(1)) .or. all(y3 == y3(1))) goto 900
  status = 3
  if (api_status /= 0 .or. .not. ieee_is_finite(r) .or. abs(r) >= 1.0_dp) goto 900

  status = 2 ! persistence failure
  call tauest_x(t2, x3, n, taux3, rhox3, persistence(1))
  call tauest_x(t2, y3, n, tauy3, rhoy3, persistence(2))
  if (any(persistence < 0)) goto 900
  call chsett4(int(block))
  call calc_t_inv_lambda
  status = 3
  if (.not. all(ieee_is_finite(t_inv_lambda))) goto 900
  call confidence(cover, selected(1))
  status = api_status
  if (status /= 0) goto 900
  status = 3
  if (.not. all(ieee_is_finite(se_r_resample2)) .or. maxval(se_r_resample2) <= 0.0_dp) goto 900
  if (.not. ieee_is_finite(r_low) .or. .not. ieee_is_finite(r_upp)) goto 900
  if (r_low > r_upp) goto 900
  corr = r
  ci(:,1) = [r_low, r_upp]
  tau = [taux3, tauy3]
  actual_block = l_mbb
  se_z = se_r_resample1
  ! Every level uses the same bootstrap curve and the same first-minimum tie rule.
  z = invtanh(r)
  do i = 2,nlevel
    one_tail = (1.0_dp - levels(i)) / 2.0_dp
    target = 1.0_dp - (2.0_dp * one_tail)
    k = minloc(abs(target - cover), dim=1)
    selected(i) = lambda(k)
    half_width = t_inv_lambda(k) * se_z
    ci(:,i) = [tanh(z + half_width), tanh(z - half_width)]
  end do
  if (.not. all(ieee_is_finite(ci))) goto 900
  degenerate = api_degenerate
  zero_inner = count(se_r_resample2 == 0.0_dp)
  z1 = r_resample1
  se2 = se_r_resample2
  grid = lambda
  status = 0

900 continue
  ! A failed multi-object allocation can leave part of the workspace allocated.
  call deallocate_resample_data
  if (allocated(t2)) deallocate(t2)
  if (allocated(x2)) deallocate(x2)
  if (allocated(y2)) deallocate(y2)
  if (allocated(x3)) deallocate(x3)
  if (allocated(y3)) deallocate(y3)
  if (allocated(x3_resample1)) deallocate(x3_resample1)
  if (allocated(y3_resample1)) deallocate(y3_resample1)
  if (allocated(x3_resample2)) deallocate(x3_resample2)
  if (allocated(y3_resample2)) deallocate(y3_resample2)
  if (allocated(rp_starts)) deallocate(rp_starts)
  nullify(api_draw)
  api_rng_context = c_null_ptr
  api_call = .false.
end subroutine p3_python
