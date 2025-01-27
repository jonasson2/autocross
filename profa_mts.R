library(MTS)
library(autocross)
var_step_series <- function(xy) {
  t = c()
  x = c()
  y = c()
  k = 1
  dk = 1
  while (k < length(xy[,1])) {
    t = c(t, k)
    x = c(x, xy[k, 1])
    y = c(y, xy[k, 2])
    k = k + dk
    dk = 10 - dk  # add 1 and 2 alternately
  }
  list(t=t, x=x, y=y)
}


# Prófum VARMA(1,1)
n = 2000
p = 1
q = 1
A = matrix(c(0.4, 0.3, 0.3, 0.4), 2, 2)
B = matrix(c(0.2, -0.5, 0.4, -0.3), 2, 2)
Sig = matrix(c(2, 1, 1, 2), 2, 2)
result = VARMAcov(Phi=A, Theta=B, Sigma=Sig, lag=0)
#th_corr=result$ccm
#print(th_corr)
result = VARMAsim(n, arlags=p, malags=q, phi=A, theta=B, sigma=Sig)
xy = result$series
vss = var_step_series(xy)
t = vss$t
x = vss$x
y = vss$y
# sim_corr = cor(x,y)
# print(paste('Simulateded correlation=', sim_corr))
print(paste('length(x)=', length(x)))
print(cor(x, y))
#result = estimate_CI(t, x, y)
#print(result)
