module mod
  save
  integer :: a = 1
end
program prog
  use mod
  print *, a
end
