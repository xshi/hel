CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C     Created by SHI Xin	  Jun 25 , 2003 		    C
C    Oct 15 , 2003 -----  Add the correlation efficient rho         C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      PROGRAM MAIN
      IMPLICIT REAL(A-H,O-Z) 


      print*,'Total sys errors:'

      sxp1=0.04; syp1=0.07  ;  rp1=1   
      sxp2=0.01; syp2=0.06  ;  rp2=1    
      sxp3=0.18; syp3=0.05  ;  rp3=0.24

      sxk1=0.04; syk1=0.08  ;  rk1=1    
      sxk2=0.10; syk2=0.14  ;  rk2=1    
      sxk3=0.13; syk3=0.08  ;  rk3=-0.24

      sxp=sqrt(sxp1**2+sxp2**2+sxp3**2)
      syp=sqrt(syp1**2+syp2**2+syp3**2)
      rhosysp=(rp1*sxp1*syp1 + rp2*sxp2*syp2 + rp3*sxp3*syp3)/(sxp*syp)

      print*,'sxp=',sxp,'   syp=',syp ,'   rhosysp=',rhosysp

      sxk=sqrt(sxk1**2+sxk2**2+sxk3**2)
      syk=sqrt(syk1**2+syk2**2+syk3**2)
      rhosysk=(rk1*sxk1*syk1 + rk2*sxk2*syk2 + rk3*sxk3*syk3)/(sxk*syk)
      print*,'sxk=',sxk,'   syk=',syk ,'   rhosysk=',rhosysk

      print*,'____________________________________________________'
      print*,'Combined errors:'

      xp=2.16; xpstat=0.64; xpsys=sxp
      yp=3.37; ypstat=1.07; ypsys=syp
               rpstat=0.96; rpsys=rhosysp

      xk=2.02; xkstat=0.54; xksys=sxk
      yk=2.83; ykstat=0.80; yksys=syk
               rkstat=0.94; rksys=rhosysk

C     combine the errors
      sigxp = sqrt(xpstat**2+xpsys**2)
      sigyp = sqrt(ypstat**2+ypsys**2)
      rhop  = (xpstat * ypstat * rpstat +
     $         xpsys * ypsys * rpsys)/(sigxp*sigyp)

      sigxk = sqrt(xkstat**2+xksys**2)
      sigyk = sqrt(ykstat**2+yksys**2)
      rhok  = (xkstat * ykstat * rkstat +
     $         xksys * yksys * rksys)/(sigxk*sigyk)

      print*,'For Pion: '
      print*,'xp=', xp ,'+/-', sigxp,'   yp=', yp, '+/-', sigyp,
     $     '   rhop=', rhop

      print*,'For Kaon:'
      print*,'xk=', xk ,'+/-', sigxk,'   yk=', yk, '+/-', sigyk,
     $     '   rhok=', rhok

      print*,'____________________________________________________'
      print*,'Final combined results >>>>>>>>>>>'


      x=com(xp, sigxp, xk, sigxk)
      y=com(yp, sigyp, yk, sigyk)

      sigx=sig(sigxp, sigxk)
      sigy=sig(sigyp, sigyk)

      rho12=rho(sigx,sigy,sigxp, sigyp, sigxk, sigyk, rhop, rhok)

      print*,'x=',x,'  +/-',sigx
      print*,'y=',y,'  +/-',sigy
      print*,'rho=',rho12
      
      END 
*********************************************************************
      real function com(x1,s1,x2,s2)
      real  x1,s1,x2,s2
      com=(x1/s1**2 + x2/s2**2)/(1/s1**2 + 1/s2**2)
      return
      END

      real function sig(s1,s2)
      real  s1,s2
      sig=sqrt(1/(1/s1**2+1/s2**2))
      return
      END

      real function rho(sigx, sigy, sigx1,sigy1, sigx2,sigy2, rho1,rho2)
      real  sigx, sigy, sigx1,sigy1, sigx2,sigy2, rho1, rho2
      rho=sigx*sigy/(sigx1*sigy1)*rho1+sigx*sigy/(sigx2*sigy2)*rho2
      return
      END

C      real function rhosys(sx,sy,sx1,sy1,sx2,sy2,sx3,sy3,r1,r2,r3)
C      real  sx, sy, sx1,sy1, sx2,sy2,sx3,sy3,r1, r2, r3
C      rhosys=(r1*sx1*sy1 + r2*sx2*sy2 + r3*sx3*sy3)/(sx*sy)
C      return
C      END

C      sigx=sig(sxp, sxk)
C      print*,'sigsys_x',sigx
C      sigy=sig(syp, syk)
C      print*,'sigsys_y',sigy
C
C      rhosys=sigx*sigy/(sxp*syp)*rhosysp+sigx*sigy/(sxk*syk)*rhosysk
C      print*,'rhosys=',rhosys
C      end
C      x1     =1.16
