C Make yourself clear! 
      SUBROUTINE TWOPN4C(LOGFLAG)
      INCLUDE 'rawdat.inc'
      INCLUDE 'trklst.inc'
      INCLUDE 'evlst.inc'
      INCLUDE 'mcmade.inc'
      INCLUDE 'pilot.inc'
      INCLUDE 'ppp.inc'
      INCLUDE 'nsquid.inc'
      INCLUDE 'xmchty.inc'
      INCLUDE '/home/shixin/com/inc/pidpsip.inc'
      INCLUDE '/home/shixin/com/inc/dofit.inc'
      INCLUDE '/home/shixin/com/inc/trkvec.inc'
      INCLUDE '/home/shixin/com/inc/npmp.inc' 
 
      INTEGER NGMMAX
      PARAMETER (NGMMAX=100)
      PARAMETER (PI=3.14159265359)
      INTEGER NGM,IDG(NGMMAX)
      real eg(NGMMAX),prbg(NGMMAX)
      real PRBEST(NGMMAX,6)
      real dc(NGMMAX),dg(NGMMAX),r12(NGMMAX)
      real r56(NGMMAX),emt(NGMMAX),tra(NGMMAX)

      real r1(NGMMAX),r2(NGMMAX),r3(NGMMAX)
      real r4(NGMMAX),r5(NGMMAX),r6(NGMMAX)
      common /rlayer/r1,r2,r3,r4,r5,r6

      INTEGER NGMX_E
      PARAMETER (NGMX_E=100)
      INTEGER NGM_esc,IDG_esc(NGMX_E)
      real eg_esc(NGMX_E)

      integer ntrkmax
      parameter (ntrkmax=100)
      integer idmdctrk(ntrkmax),idindx(ntrkmax)
      integer nmdctrk
      common /mdctrk/nmdctrk,idmdctrk,idindx
      
      COMMON/PAWC/Hmem(1000000)
      COMMON/USEARG/USEARG(10)
      INTEGER USEARG,LOGFLAG
      LOGICAL   IDMU1,IDMU2
      INTEGER newmuid1,newmuid2
      INTEGER IDPTID1(10)
      REAL  CHISQIDPI(4),CHISQIDKA(4),CHISQIDPR(4)
      INTEGER NDOFID(4),NDOFALL
      
      PARAMETER(NTPID = 1,NTAG = 100)
      REAL*4 TNTP(ntag)
      INTEGER ITNTP(ntag)
      EQUIVALENCE (ITNTP(1),TNTP(1))
      
      real probpion(4),probkaon(4),probpr(4)
     
      REAL SP(4),PB(4),PF(4)  ! CHANGED BY SHIXIN       
      integer isys,i1,i2,i3,iv,iflag
      real cthtv,cthtx,ctheta,phipi   !call calangle.f

      LOGFLAG=0
      call getecm(irun,ecm,energy)
      ivs = 103
      
      DO I=1,NTAG
        TNTP(I)=999.
      ENDDO
      TNTP(33)=0.
      TNTP(43)=0.
      TNTP(53)=0.
      
      DO I=1,10
        IDPTID1(I) = 0
      ENDDO

      if(irun.gt.10000.and.irun.lt.20000)then     !J/psi
        call jpsi_run_status(irun,ijk1)
        ijk2 = 0
      else
        call skip_bad_run(ijk1,ijk2)
      endif
      TNTP(13) = float(ijk1)

c------------------------------------------EVENT general info.
      
      TNTP( 1) = FLOAT(IRUN)
      TNTP( 2) = FLOAT(IREC)
      TNTP( 3) = ENERGY

c------------------------------------------photon ID
      CALL phtid_bes1(ngm,idg,prbg,prbest,eg,dc,dg,r12,r56,emt,tra)
      CALL phtid_esc(ngm_esc,idg_esc,eg_esc)
      TNTP( 4) = FLOAT(ngm_esc*100+NGM)
      
      if(NGM.LT.1)goto 800      !at least two photons
      NP13=NP13+1               !NGM.GE.1

      idggp=0
      prbgmax=0.
      do i=1,ngm
        if(eg(i).gt.prbgmax.and.eg(i).gt.-0.050)then
          prbgmax=eg(i)
          idggp=idg(i)
          indxgp=i
        endif
      enddo
      if(idggp.eq.0)goto 800
      
C-----------------------------------------select TWO good charged tracks
C The tracks are arranged in the following way: + -
C
      IK    = 0
      DO 60 I=1,NMDCTRK
        IF(PP(I).LT.0.0.OR.PP(I).GT.3.0)GOTO 60
        IK=IK+1
        if(ich(I).eq.+1)THEN
          idpip=i
        else
          idpim=i
        endif
   60 CONTINUE
      IF(IK.NE.2)GOTO 800       !Two good charged tracks needed
      NP14=NP14+1               !IK.EQ.2
      
C---------------------------------------------PI PI L L
      CALL MYMUID(idmdctrk(IDPIP),NEWMUID1,IDMU1,NEWMUID2,IDMU2)
      TNTP(5) = NEWMUID2
      CALL MYMUID(idmdctrk(IDPIM),NEWMUID1,IDMU1,NEWMUID2,IDMU2)
      TNTP(6) = NEWMUID2
C
      if(TNTP(5)+TNTP(6).ge.5)goto 800      !dimu event
      NP15=NP15+1               !NO dimu MUIDP+MUIDM<5
C
      EG1=0.
      ISUB5  = TRKSUB(5,idmdctrk(IDPIP))
      IF(ISUB5.NE.0)EG1=TRK(ISUB5+14)
      TNTP(7) = EG1
      EG1=0.
      ISUB5  = TRKSUB(5,idmdctrk(IDPIM))
      IF(ISUB5.NE.0)EG1=TRK(ISUB5+14)
      TNTP(8) = EG1
C
C      if(TNTP(7).gt.1.0.and.TNTP(8).gt.1.0)goto 800      !bhabha event
C      if(TNTP(7).gt.1.0.or.TNTP(8).gt.1.0)goto 800      ! tight cut bhabha event
      NP16=NP16+1               !NO bhabha     
C
Csx------------------------------------------------
C check dEdX 
      dedx=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIP))
      IF(ISUB7.NE.0)dedx=TRK(ISUB7+11)! run correction for pulse heights
      TNTP(46) = dedx ! dedxp 
      dedx=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIM))
      IF(ISUB7.NE.0)dedx=TRK(ISUB7+11) ! run correction for pulse heights
      TNTP(47) = dedx ! dedxm
      
      dedx=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIP))
      IF(ISUB7.NE.0)dedx=TRK(ISUB7+8)! pulse height from truncated mean
      TNTP(48) = dedx ! rdedxp 
      dedx=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIM))
      IF(ISUB7.NE.0)dedx=TRK(ISUB7+8) 
      TNTP(49) = dedx ! rdedxm
Csx------------------------------------------------
C use XSE to cut e+e- background 
      XSE1=0.
      XSMU1=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIP))
      IF(ISUB7.NE.0)XSE1=TRK(ISUB7+15) 
      TNTP(14) = XSE1 !  of sigma from e
      IF(ISUB7.NE.0)XSMU1=TRK(ISUB7+16) 
      TNTP(15) = XSMU1

      XSE2=0.
      XSMU2=0.
      ISUB7  = TRKSUB(7,idmdctrk(IDPIM))
      IF(ISUB7.NE.0)XSE2=TRK(ISUB7+15) 
      TNTP(19) = XSE2 !  of sigma from e
      IF(ISUB7.NE.0) XSMU2=TRK(ISUB7+16) 
      TNTP(20) = XSMU2

Csx------------------------------------------------

C---------------------------------------------DO Partical ID
      npion=0
      nkaon=0
      nprot=0
      DO ITKPAR=1,nmdctrk
        ndofid(itkpar)=0
        chisqidpi(itkpar)=0.
        chisqidka(itkpar)=0.
        chisqidpr(itkpar)=0.
      ENDDO
      DO I=1,NMDCTRK
         ITKPAR = idmdctrk(i)
C>>>>>>> variant B
         call tofdxpid
         probpion(i)=PROBab(3)
         probkaon(i)=PROBab(4)
         probpr(i)  =PROBab(5)
         if(itfgood.eq.1)then
            
            ndofid(i)=ndofid(i)+1
            chisqidpi(i)=chisqidpi(i)+tofchi(3)**2
            chisqidka(i)=chisqidka(i)+tofchi(4)**2
            chisqidpr(i)=chisqidpr(i)+tofchi(5)**2
C     sx------------------------------------------------
C     check TOF 
C     TNTP(36)=tofchi(3) ! tofchipi
C     TNTP(37)=tofchi(4) ! tofchika
C     TNTP(38)=tofchi(5) ! tofchipr
C     TNTP(39)=trk(trksub(4,ITKPAR)+7) !pulse
C     sx------------------------------------------------
         endif
         if(idxgood.eq.1)then
            ndofid(i)=ndofid(i)+1
            chisqidpi(i)=chisqidpi(i)+dedxchi(3)**2
            chisqidka(i)=chisqidka(i)+dedxchi(4)**2
            chisqidpr(i)=chisqidpr(i)+dedxchi(5)**2
         endif
         IF(idpart.eq.3.and.probpion(i).gt.0.01)THEN
            
            IDPTID1(I)=3        !Pion
            npion = npion + 1
         ENDIF
         IF(idpart.eq.4.and.probkaon(i).gt.0.01)THEN
            IDPTID1(I)=4        !kaon
            nkaon = nkaon + 1
         ENDIF                                   
         if(idpart.eq.5.and.probpr(i).gt.0.01)then
            IDPTID1(I)=5        !proton
            nprot = nprot + 1
C======= end
        endif
      ENDDO

C>>>>>>> variant B
C cosmic ray veto |TOF1-TOF2|<4ns 
         ITKPAR1 = idmdctrk(1)
         ITKPAR2 = idmdctrk(2)
         
         TNTP(36)=trk(trksub(4,ITKPAR1)+4)! tof1
         TNTP(37)=trk(trksub(4,ITKPAR2)+4)! tof2
      
C======= end
C-- redefine id** as trklst index
      idpip1 = idmdctrk(idpip)
      idpim1 = idmdctrk(idpim)
      idxpip1 = idindx(idpip1)
      idxpim1 = idindx(idpim1)
      
      ndofall=4+ndofid(1)+ndofid(2)
      TNTP(9) = ndofall
      
      eescmax=-999.
      do i=1,ngm_esc
         if(eescmax.lt.eg_esc(i))then
            eescmax=eg_esc(i)
         endif
      enddo
      TNTP(10) = eescmax
      
      TNTP(11) = eg(indxgp)
      TNTP(12) = prbg(indxgp)

C--------------------------------------------fill event vertex
      NVERTEX= 1                ! CUT FOR VERTEX
      CALL VXKWIK(ISTV,NVERTEX)
      IVERT2 = VRTSUB(2, 1)
      IF(ISTV.GE.2.OR.IVERT2.LE.0)THEN
         XIR  = 9.0
         YIR  = 9.0
         ZIR  = 9.0
      ELSE
         XIR  = TRK(IVERT2 + 5)
         YIR  = TRK(IVERT2 + 6)
         ZIR  = TRK(IVERT2 + 7)
      ENDIF
      TNTP(16)=XIR
      TNTP(17)=YIR
      TNTP(18)=ZIR
      
      TNTP(21) = pp(IDxPIp1)
      TNTP(22) = probpion(IDxPIp1)
      TNTP(23) = probkaon(IDxPIp1)
      TNTP(24) = probpr(IDxPIp1)
      TNTP(25) = pz(IDxPIp1)

      TNTP(26) = pp(IDxPIm1)
      TNTP(27) = probpion(IDxPIm1)
      TNTP(28) = probkaon(IDxPIm1)
      TNTP(29) = probpr(IDxPIm1)
      TNTP(30) = pz(IDxPIM1)

C<<<<<<< variant A
C_ Add more information to calculat acol

      TNTP(56) = px(IDxPIp1)
      TNTP(57) = py(IDxPIp1)

      TNTP(58) = px(IDxPIm1)
      TNTP(59) = py(IDxPIM1)

C>>>>>>> variant B
C======= end
 210  continue
C---------------------------------------------DO 4C-fit TO PI+PI- GAMMA
      Call psinit(energy,npart,nres,rc)
      Call psdotrkc(1,xmpich,1,1,1,idpip1,npart,mass,pspec,rc,ivs)
      Call psdotrkc(2,xmpich,1,1,1,idpim1,npart,mass,pspec,rc,ivs)
      Call pstrkcp(3,0.,1,1,1,idggp,npart,mass,pspec,rc,ivs)
      FITPPG=newsquid(npart,mass,pspec,nres,
     &  mres,pmem,mem,rc,chisqPP,pfit)
      
      IF(.NOT.FITPPG.OR.CHISQPP.GT.20)GOTO 220
      NP17=NP17+1               !Fitppg=true Chisqpp<20
      Call squidpul(pspec,rc,pull)
      call cppfit2p(pfit,mass,3,p)
      
      CALL ADD$(1,2,8)
      
      TNTP(31) = chisqidpi(1)+chisqidpi(2) !CSQIDPI
      TNTP(32) = CHISQPP        !CSQ4CPI
      TNTP(33) = PROB(chisqPP,4) !PRB4CPI
      TNTP(34) = p(5,8)
      TNTP(35) = p(4,3)
C-------------------Lorentz Transformations for Pi+
      SP(1)= p(1,8)
      SP(2)= p(2,8)
      SP(3)= p(3,8)
      SP(4)= p(4,8)
      SM   = p(5,8)
     
      PB(1)= p(1,1)
      PB(2)= p(2,1)
      PB(3)= p(3,1)
      PB(4)= p(4,1) 

      i1 = 1
      i2 = 2
      i3 = 0
      iv = 3
      iflag = 0
      isys =  9
      do i=1,4
        if(i.ge.3)then
          p(i,isys)=energy/2.
        else
          p(i,isys)=0.0
        endif
      enddo

      call angle(isys,i1,i2,i3,iv,cthtv,cthtx,ctheta,phipi,iflag)

      thetaG=ACOS(p(3,3)/SQRT(p(1,3)**2+p(2,3)**2+p(3,3)**2))

      TNTP(91)=thetaG           !thetaGP
      TNTP(92)=ACOS(ctheta)     !thetaP
      TNTP(93)=phipi            !phiP


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
       LOGFLAG = LOGFLAG + 1
       
 
C---------------------------------------------DO 4C-fit to K K GAMMA
 220  continue   

      Call psinit(energy,npart,nres,rc)
      Call psdotrkc(1,xmkch,1,1,1,idpip1,npart,mass,pspec,rc,ivs)
      Call psdotrkc(2,xmkch,1,1,1,idpim1,npart,mass,pspec,rc,ivs)
      Call pstrkcp(3,0.,1,1,1,idggp,npart,mass,pspec,rc,ivs)
      FITKKG=newsquid(npart,mass,pspec,nres,
     &  mres,pmem,mem,rc,chisqKK,pfit)

      IF(.NOT.FITKKG.OR.CHISQKK.GT.20)GOTO 250
      NP18=NP18+1     
      Call squidpul(pspec,rc,pull)
      call cppfit2p(pfit,mass,3,p)
      
      CALL ADD$(1,2,8)

      TNTP(41) = chisqidka(1)+chisqidka(2)
      TNTP(42) = CHISQKK
      TNTP(43) = PROB(chisqKK,4)
      TNTP(44) = p(5,8)
      TNTP(45) = p(4,3)

C-------------------Lorentz Transformations for K+
    
      SP(1)= p(1,8)
      SP(2)= p(2,8)
      SP(3)= p(3,8)
      SP(4)= p(4,8)
      SM   = p(5,8)
     
      PB(1)= p(1,1)
      PB(2)= p(2,1)
      PB(3)= p(3,1)
      PB(4)= p(4,1) 

      i1 = 1
      i2 = 2
      i3 = 0
      iv = 3
      iflag = 0
      isys =  9
      do i=1,4
        if(i.ge.3)then
          p(i,isys)=energy/2.
        else
          p(i,isys)=0.0
        endif
      enddo

      call angle(isys,i1,i2,i3,iv,cthtv,cthtx,cthetak,phik,iflag)
      thetaG=ACOS(p(3,3)/SQRT(p(1,3)**2+p(2,3)**2+p(3,3)**2))

      TNTP(88)=thetaG           !thetaGK
      TNTP(89)=ACOS(cthetak)     !thetaK
      TNTP(90)=phik            !phiK

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      LOGFLAG = LOGFLAG + 1
      
C---------------------------------------------DO 4C-fit to P P gamma
 250  continue

CC      Call psinit(energy,npart,nres,rc)
CC      Call psdotrkc(1,xmp,1,1,1,idpip1,npart,mass,pspec,rc,ivs)
CC      Call psdotrkc(2,xmp,1,1,1,idpim1,npart,mass,pspec,rc,ivs)
CC      Call pstrkcp(3,0.,1,1,1,idggp,npart,mass,pspec,rc,ivs)
CC      FITPRG=newsquid(npart,mass,pspec,nres,
CC     &  mres,pmem,mem,rc,chisqPR,pfit)
CC      
CC      IF(.NOT.FITPRG.OR.CHISQPR.GT.20)GOTO 299
CC      NP19=NP19+1     
CC      Call squidpul(pspec,rc,pull)
CC      call cppfit2p(pfit,mass,3,p)
CC
CC      CALL ADD$(1,2,8)
CC
CC      TNTP(51) = chisqidpr(1)+chisqidpr(2)
CC      TNTP(52) = CHISQPR
CC      TNTP(53) = PROB(chisqPR,4)
CC      TNTP(54) = p(5,8)
CC      TNTP(55) = p(4,3)
CC
CC      LOGFLAG = LOGFLAG + 1
CC      
CC 299  continue

C............................................................Save events     
      
      if(LOGFLAG.ne.0)then
        CALL HFN(NTPID, TNTP)
      endif
      
 800  RETURN
      END

C********************************************************************
      SUBROUTINE BESEVN(IIREC,EVNARG)

      INCLUDE 'rawdat.inc'
      INCLUDE 'trklst.inc'
      INCLUDE 'evlst.inc'
      INCLUDE 'mcmade.inc'
      INCLUDE '/home/shixin/com/inc/trkvec.inc'
      INCLUDE '/home/shixin/com/inc/npmp.inc'

      integer ntrkmax
      parameter (ntrkmax=100)
      integer idmdctrk(ntrkmax),idindx(ntrkmax)
      integer nmdctrk
      common /mdctrk/nmdctrk,idmdctrk,idindx
      
      COMMON/PAWC/Hmem(1000000)
      INTEGER*4 LOGFLAG
      COMMON/USEARG/USEARG(10)
      INTEGER USEARG
      LOGFLAG=0

      IF(IREC/10000.EQ.FLOAT(IREC)/10000.)
     &  WRITE(6,*)'Processing RUN, REC',IRUN,IREC
CYUANCZ ----------------------------------------------- SKIP BAD RUN
      NP0=NP0+1
      if(irun.gt.10000.and.irun.lt.20000)then     !J/psi
        call jpsi_run_status(irun,istatus)
        CALL HF1(2,float(istatus),1.)
c        if(istatus.ne.1)goto 800
      else
        call skip_bad_run(ijk1,ijk2)
        CALL HF1(2,float(ijk1),1.)
c        if(ijk1.le.1)goto 800
      endif
      NP1=NP1+1

      CALL HF1(3,FLOAT(NCHRG),1.)
      NMDCTRK = 0
      do i=1,ntrkmax
        IDMDCTRK(i) = 0
        IDINDX(i) = 0
      enddo
      do i=1,nchrg
        isub1=trksub(1,i)
        if(isub1.gt.0)then
          ITYPE = ITRK(ISUB1+5)
          IDC = IBITS(ITYPE,0,1)
          IF(IDC.NE.0)THEN      !MDC-track
            NMDCTRK = NMDCTRK + 1
            idmdctrk(NMDCTRK) = i
            idindx(i)=NMDCTRK
          ENDIF
        endif
      enddo
      IF(NMDCTRK.NE.2) GOTO 800 !exactly 2 MDC tracks
      NP2=NP2+1
      CALL HF1(4,FLOAT(NNEU),1.)
      IF(NNEU.LT.0)GOTO 800     !nneu
      NP3=NP3+1
      
      DO I=1,NMDCTRK            !charged tracks
        ISUB3= TRKSUB(3,IDMDCTRK(i))
        IF(ISUB3.LE.0)GOTO 800
        IF(IABS(ITRK(ISUB3+2)).NE.1)GOTO 800 !ICHRG
      ENDDO
      NP4=NP4+1

      DO I=1,NMDCTRK
         ISUB3= TRKSUB(3,IDMDCTRK(i))
        CALL HF1(5,FLOAT(ITRK(ISUB3+14)),1.)
        IF(ITRK(ISUB3+14).NE.2) GOTO 800 !MFIT
C     &       .AND.ITRK(ISUB3+14).NE.-19
C     &       .AND.ITRK(ISUB3+14).NE.-9)GOTO 800 !MFIT
        
      ENDDO
      NP5=NP5+1
      
      ITOTCH= 0
      DO I=1,NMDCTRK
         ISUB3= TRKSUB(3,IDMDCTRK(i))
         ITOTCH = ITOTCH + ITRK(ISUB3+2)
      ENDDO
      CALL HF1(6,FLOAT(ITOTCH),1.)
      IF(ITOTCH.NE.0)GOTO 800   !TOT CHRG
      NP6=NP6+1
      
      NP7=NP7+1
C................................read vectors from TRKLST
      DO 10 I=1,NMDCTRK
         ISUB3= TRKSUB(3,IDMDCTRK(i)) !TRK FIT,MDC
         ICH(I)=ITRK(ISUB3+2)
         PXY(I)=TRK(ISUB3+3)
         PX(I)=TRK(ISUB3+4)
         PY(I)=TRK(ISUB3+5)
         PZ(I)=TRK(ISUB3+6)
         PP(I)=TRK(ISUB3+7)
         MFIT(I)=ITRK(ISUB3+14)
 10   CONTINUE
      
      CALL TWOPN4C(LOGFLAG)
      CALL HF1(7,FLOAT(LOGFLAG),1.)
      IF(LOGFLAG.EQ.0)GOTO 800
      NP8=NP8+1
      
      CALL EVLOG(II)

 800  RETURN
      END
      

C************************************************************************   
      SUBROUTINE BESIOR(IRUN,IORARG)
      INCLUDE 'pilot.inc'
      INCLUDE '/home/shixin/com/inc/npmp.inc'
      COMMON/USEARG/USEARG(10)
      INTEGER USEARG
      INTEGER IORARG(10)
C      COMMON/ICNT/ICNT(10)
C      INTEGER ICNT(10)
C      DATA ICNT /10*0/
      COMMON/PAWC/Hmem(1000000)
     
      DATA NP0,NP1,NP2,NP3,NP4/5*0/
      DATA NP5,NP6,NP7,NP8,NP9/5*0/
      
      DATA NP10,NP11,NP12,NP13,NP14/5*0/
      DATA NP15,NP16,NP17,NP18,NP19/5*0/
      
      PARAMETER(NTPID = 1, NTAG =100)
      CHARACTER TAGS(NTAG)*8
      DATA NPRIME/1000/
      DATA TAGS  / '    IRUN','    IREC','     ECM',' NESCNGM', 
     &     '   MUIDP','   MUIDM','    SCEP','    SCEM','  NDFALL',
     &     ' EESCMAX',
     &     '      EG','    PRBG','    IJK1','    XSE1','   XSMU1',
     &     '     XIR','     YIR','     ZIR','    XSE2','   XSMU2',!20
C partid PIP
     &     '     PPP','  PRBPIP','  PRBKAP','  PRBPRP','     PZP',
C     partid PIM
     &  '     PPM','  PRBPIM','  PRBKAM','  PRBPRM','     PZM',   !30
C 4C-fit to PiP PiM gamma 
     &     ' CSQIDPI',' CSQ4CPI',' PRB4CPI','  MALLPI','    EGPI',
C<<<<<<< variant A
C     &     'tofchipi','tofchika','tofchipr','   pulse','   ICNT5',!40
C>>>>>>> variant B
     &     '    tof1','    tof2','tofchipr','   pulse','   ICNT5',!40
C======= end
C     4C-fit to KaP KaM gamma 
     &     ' CSQIDKA',' CSQ4CKA',' PRB4CKA','  MALLKA','    EGKA',
     &     '   dedxp','   dedxm','  rdedxp','  rdedxm','  ICNT10',!50
C     4C-fit to PrP PrM gamma 
     &     ' CSQIDPR',' CSQ4CPR',' PRB4CPR','  MALLPR','    EGPR',
C<<<<<<< variant A
C    More info about the two charges
     &  '     PXP','     PYP','     PXM','     PYM','     T60',   !60
C>>>>>>> variant B
C     &  '     T56','     T57','     T58','     T59','     T60',   !60
C======= end
      
C 61-65
     &     ' P4CPIPX',' P4CPIPY',' P4CPIPZ','  E4CPIP','  M4CPIP',
C     66-70
     &     '     PFX','     PFY','     PFZ','     PFE','     T70',
C     71-75
     &     '   P4CGX','   P4CGY','   P4CGZ','    E4CG','    M4CG',
C     76-80
     &     ' P4CPIMX',' P4CPIMY',' P4CPIMZ','  E4CPIM','  M4CPIM',
C     81-85
     &     ' PALLPIX',' PALLPIY',' PALLPIZ','  EALLPI',' MALLPI2',
C     86-90
     &     '     thP','     phP',' thetaGK','  thetaK','    phiK',
C     91-95
     &     ' thetaGP','  thetaP','    phiP','     T95','     T96',
C     96-100
     &     '     T96','     T97','     T98','     T99','    T100'/
      



      DO IY=1,10
         USEARG(IY) = IORARG(IY)
      ENDDO
      
      CALL HLIMIT(1000000)
      OPEN(80,FORM='UNFORMATTED',RECL=4096,
     &     ACCESS='DIRECT',STATUS='NEW')
      CALL HRFILE(80,'TWOPN','N')
*     Ntuple:
      CALL HBOOKN(NTPID,'TWOPN ntuple',NTAG,
     &     'TWOPN',NPRIME,TAGS)
      CALL HBOOK1(2,'IJK1',20,-10.0,9.5,0.)
      CALL HBOOK1(3,'NCHRG',20,-0.5,19.5,0.)
      CALL HBOOK1(4,'NNEU',20,-0.5,19.5,0.)
      CALL HBOOK1(5,'MFIT',40,-19.5,20.5,0.)
      CALL HBOOK1(6,'Net Charge',5,-2.5,2.5,0.)
      CALL HBOOK1(7,'Log Flag',4,-0.5,3.5,0.)
      
      RETURN
      END

C********************************************************************
      SUBROUTINE BESEOR(IRUN,EORARG)
      COMMON/PAWC/Hmem(1000000)
      PARAMETER (NTPID = 1,NTAG = 100)
      INCLUDE '/home/shixin/com/inc/npmp.inc'
      
      WRITE(6,*)'---------------------------------------------------'
      WRITE(6,*)'Total events:         ',NP0
      WRITE(6,*)'1. NMDCTRK=2:         ',NP2
      WRITE(6,*)'2. ISUB3>0,ICHRG=+-1: ',NP4
      WRITE(6,*)'3. MFIT=2:            ',NP5
      WRITE(6,*)'4. Q1+Q2=0:           ',NP6
      WRITE(6,*)'5. NGM>=1:            ',NP13
      WRITE(6,*)'6. P<=3.0GeV:         ',NP14
      WRITE(6,*)'7. MUIDP+MUIDM<5:     ',NP15
      WRITE(6,*)'8. ESCP OR ESCM<=1.0: ',NP16
      WRITE(6,*)'9. Pass Kinematic Fit:',NP8
      WRITE(6,*)'---------------------------------------------------'
      WRITE(6,*)'---------------------------------------------------'
      WRITE(6,*)'IJK .GT. 0  :      ',NP1
      WRITE(6,*)'NNEU>=0:           ',NP3
      WRITE(6,*)'CHISQPP<20:pi      ',NP17
      WRITE(6,*)'IFIRLAY.ne.6:      ',NP10
      WRITE(6,*)'CSEMIT>0.8  :      ',NP11
      WRITE(6,*)'ACOSD>15:          ',NP12
      WRITE(6,*)'K:                 ',NP18
      WRITE(6,*)'p:                 ',NP19
      WRITE(6,*)'---------------------------------------------------'

      CALL HROUT(0,ICYCLE,'T')
      CALL HREND('TWOPN')

      RETURN
      END
