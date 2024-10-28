module test
contains
  subroutine sub(s)
    character*(*), intent(in) :: s
    print '(A)', '<'//s//'>'
  end subroutine sub
end module test
program prufa
  use test
  call sub('abc')
  print '(lzp,"<"F0.2">")', 0.14159
end program prufa
