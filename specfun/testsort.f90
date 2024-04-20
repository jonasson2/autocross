program main
  use sort_module
  implicit none
  integer, parameter :: N = 10
  double precision :: arr(N) = [5.0, 1.0, 3.0, 8.0, 2.0, 7.0, 4.0, 9.0, 6.0, 10.0]
  integer :: i
  call sort(arr)
  print *, "Sorted array:"
  do i = 1, N
    print *, arr(i)
  end do
end program main
