  subroutine ss(x, s)
    real(8), intent(in) :: x(2)   ! 1D array of the first two elements as input
    real(8), intent(out) :: s     ! Output variable for the sum
    s = x(1) + x(2)
    print *,'s=', s
end subroutine ss
