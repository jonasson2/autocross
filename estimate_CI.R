function estimate_CI(data, confidence_level) result(ci)
  implicit none
  real, intent(in) :: data(:)
  real, intent(in) :: confidence_level
  real :: ci(2)
  
  integer :: n
  real :: mean, std_dev, t_value
  
  ! Calculate sample size
  n = size(data)
  if (n <= 1) then
  print*, "Error: Not enough data to estimate confidence interval"
  return
  end if
  
  ! Calculate mean
  mean = sum(data) / n
  
  ! Calculate standard deviation
  std_dev = sqrt(sum((data - mean)**2) / (n - 1))
  
  ! Calculate critical t-value (for simplicity, using a pre-defined value or approximate formula)
  ! Note: In practice, you would use a t-distribution table or an external library
  t_value = 1.96  ! Approximate for 95% confidence level and large sample size
  
  ! Calculate confidence interval
  ci(1) = mean - t_value * (std_dev / sqrt(n))
  ci(2) = mean + t_value * (std_dev / sqrt(n))

end function estimate_CI