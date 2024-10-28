      SUBROUTINE avevar(data,n,ave,var)
      INTEGER n
      DOUBLE PRECISION ave,var,data(n)
      INTEGER j
      DOUBLE PRECISION s,ep
      ave=0.0d0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0d0
      ep=0.0d0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      return
      END
