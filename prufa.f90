subroutine prufa(n,x)
  integer, intent(in) :: n
  double precision, intent(out) :: x(2,3)
  print *, 'n=', n
  x = reshape([1, 2, 3, 4, 5, 6], shape = [2,3])
end subroutine prufa
