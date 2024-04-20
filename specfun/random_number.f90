program test_random
    implicit none
    double precision :: x
    integer :: i

    ! Seed the random number generator
    call random_seed()

    ! Generate and print 10 standard normal random numbers
    do i = 1, 10
        call random_number(x)
        print *, 'Random number:', x
    end do
end program test_random
