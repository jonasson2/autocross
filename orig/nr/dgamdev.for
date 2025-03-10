      FUNCTION gamdev(ia,idum)
      INTEGER ia,idum
      DOUBLE PRECISION gamdev
CU    USES ran1
      INTEGER j
      DOUBLE PRECISION am,e,s,v1,v2,x,y,ran1
      if(ia.lt.1)pause 'bad argument in gamdev'
      if(ia.lt.6)then
        x=1.d0
        do 11 j=1,ia
          x=x*ran1(idum)
11      continue
        x=-log(x)
      else
1         v1=ran1(idum)
          v2=2.d0*ran1(idum)-1.d0
        if(v1**2+v2**2.gt.1.d0)goto 1
          y=v2/v1
          am=ia-1
          s=sqrt(2.d0*am+1.d0)
          x=s*y+am
        if(x.le.0.d0)goto 1
          e=(1.d0+y**2)*exp(am*log(x/am)-s*y)
        if(ran1(idum).gt.e)goto 1
      endif
      gamdev=x
      return
      END
