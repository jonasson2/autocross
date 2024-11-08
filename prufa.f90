subroutine prufa(n, x, y, sum)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n)
  double precision, intent(out) :: sum
  print *, 'n=', n
  sum = x(1) + y(1)
  return
end subroutine prufa
