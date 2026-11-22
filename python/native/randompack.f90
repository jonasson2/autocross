! The Python wheel exposes integer draws through Rng.int, not a Fortran
! library. These bindings satisfy the legacy branches of p3_modules; the
! Python entry point uses its draw callback exclusively.
module randompack
  use iso_c_binding, only: c_int32_t
  implicit none
  private
  public :: randompack_rng

  type :: randompack_rng
  contains
    procedure, nopass :: create => legacy_create
    procedure, nopass :: seed => legacy_seed
    procedure, nopass :: int => legacy_int
  end type
contains
  subroutine legacy_create(engine)
    character(len=*), intent(in) :: engine
    error stop 'Legacy RNG creation is not available in the Python extension'
  end subroutine

  subroutine legacy_seed(seed)
    integer(c_int32_t), intent(in) :: seed
    error stop 'Legacy RNG seeding is not available in the Python extension'
  end subroutine

  subroutine legacy_int(values, lower, upper)
    integer(c_int32_t), intent(out) :: values(:)
    integer, intent(in) :: lower, upper
    error stop 'Legacy RNG draws are not available in the Python extension'
  end subroutine
end module
