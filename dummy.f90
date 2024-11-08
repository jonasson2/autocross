subroutine pearsont3sub(n, t1, x1, y1, alpha, corr, ci, taux, tauy)
  integer, parameter :: dp = kind(1d0)
  integer, intent(in) :: n
  real(dp), intent(in) :: t1(n), x1(n), y1(n), alpha
  real(dp), intent(out) :: corr, ci(2), taux, tauy
  print *, 'n=', n
  print *, 't1=', t1
  corr = 0.5
  ci = [0.4, 0.6]
  taux = 3
  tauy = 4
  return
end subroutine pearsont3sub
