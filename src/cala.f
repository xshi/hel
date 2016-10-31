      PROGRAM MAIN
      IMPLICIT REAL(A-H,O-Z) 



      x=2.08; sx=0.44
      y=3.03; sy=0.66

      x1=x-sx; x2=x+sx
      y1=y-sy; y2=y+sy
      
      A=(1+y**2-2*x**2)/(1+y**2+2*x**2)
      A1=(1+y1**2-2*x1**2)/(1+y1**2+2*x1**2)
      A2=(1+y2**2-2*x2**2)/(1+y2**2+2*x2**2)
      print*,'A=',A
      print*,'A1=',A1
      print*,'A2=',A2
      
      END 
