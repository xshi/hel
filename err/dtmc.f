      PROGRAM MAIN
      IMPLICIT REAL(A-H,O-Z) 

C _{pi}
      sxp1=0.14 ;  sxp2=0.12 
      syp1=0.0  ;  syp2=0.05
       rp1=1.0  ;   rp2=0.37
c _{K}
      sxk1=0.04  ;  sxk2=0.12 
      syk1=-0.06 ;  syk2=0.05
       rk1=1.0   ;   rk2=0.01

c      sxp1=0.12 ;  sxp2=0.12 
C      syp1=0.0  ;  syp2=0.05
C       rp1=1.0  ;   rp2=0.34
Cc _{K}
C      sxk1=0.01  ;  sxk2=0.14 
C      syk1=-0.04 ;  syk2=0.05
C       rk1=1.0   ;   rk2=0.13
C


      sxp=sqrt(sxp1**2+sxp2**2)
      syp=sqrt(syp1**2+syp2**2)
      rhop=(sxp1*syp1*rp1+sxp2*syp2*rp2)/(sxp*syp)

      sxk=sqrt(sxk1**2+sxk2**2)
      syk=sqrt(syk1**2+syk2**2)
      rhok=(sxk1*syk1*rk1+sxk2*syk2*rk2)/(sxk*syk)

      print*,'sxp=', sxp 
      print*,'syp=', syp 
      print*,'rhop=', rhop

      print*,'KKK'
      print*,'sxk=', sxk 
      print*,'syk=', syk 
      print*,'rhok=', rhok

       END 
