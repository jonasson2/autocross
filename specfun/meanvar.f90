! Use Welford's algorithm, cf. Wikipedia
module meanvar_module
contains
  subroutine meanvar(x, mean, var)
    implicit none
    double precision, dimension(:), intent(in):: x
    double precision, intent(out):: mean, var
    double precision delta, M2
    integer :: i, n
    n = size(x)
    mean = 0.0
    M2 = 0.0
    do i = 1, n
      delta = x(i) - mean
      mean = mean + delta/i
      M2 = M2 + delta*(x(i) - mean)
    end do
    if (n > 1) then
      var = M2 / (n - 1)
    else
      var = 0.0
    endif
  end subroutine meanvar
end module meanvar_module
