subroutine prufa(n, x, y, sum)
  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n)
  double precision, intent(out) :: sum
  print *, 'n=', n
  print *, 'x,y=', shape(x), shape(y)
  sum = 2 + 2
  return
end subroutine prufa
