clear y
t=linspace(1,10);
n = 5;
for a=1:5
  y(:,a) = gammainc(t,a);
end
plot(y, linew=2)
title('Regularized incomplete gamma function')
xlabel('$$\frac{1}{\gamma(a)}\int_0^x t^{a-1}e^{-t}dt$$', ...
  'interpreter', 'latex', FontSize=18)
grid
legend("a=" + num2str((1:5)'))