      SUBROUTINE DISENT(NEV,S,NFL,USER,CUTS,USD1,USD2,ORDER,SCALE)
      IMPLICIT NONE
      INTEGER NEV, NFL, USD1, USD2, ORDER
      DOUBLE PRECISION S, CF, CA, TR, SCALE, NPOW1, NPOW2, CUTOFF, SCL
      EXTERNAL USER, CUTS
      LOGICAL SCALE_VAR

      CF = 4d0/3d0
      CA = 3d0
      TR = 0.5d0
      NPOW1 = 2D0
      NPOW2 = 4D0
      CUTOFF = 1d-8
      SCALE_VAR = .false.

      CALL DISENTFULL(NEV,S,NFL,USER,CUTS,USD1,USD2,NPOW1,NPOW2,CUTOFF
     $     ,SCALE,ORDER,CF,CA,TR,SCALE_VAR)
      END

      SUBROUTINE DISENTEXTENDED(NEV,S,NFL,USER,CUTS,USD1,USD2,ORDER
     $     ,SCALE,CF,CA,TR,SCALE_VAR)
      IMPLICIT NONE
      INTEGER NEV, NFL, USD1, USD2, ORDER
      DOUBLE PRECISION S, CF, CA, TR, SCALE, NPOW1, NPOW2, CUTOFF
      EXTERNAL USER, CUTS
      LOGICAL SCALE_VAR

      NPOW1 = 2D0
      NPOW2 = 4D0
      CUTOFF = 1d-8
!      SCALE_VAR = .false.

      CALL DISENTFULL(NEV,S,NFL,USER,CUTS,USD1,USD2,NPOW1,NPOW2,CUTOFF
     $     ,SCALE,ORDER,CF,CA,TR,SCALE_VAR)
      END
      
      SUBROUTINE DISENTFULL(NEV,S, NFL,USER,CUTS,USD1,USD2,NPOW1,NPOW2
     $     ,CUTOFF_IN, SCALE_in,ORDER,CF_in,CA_in,TR_in,SCALE_VAR_in)
      IMPLICIT NONE
C---CALCULATE DIS EVENT FEATURES TO NEXT-TO-LEADING ORDER
C   ACCORDING TO THE METHOD OF CATANI AND SEYMOUR NPB485 (1997) 291
C
C   VERSION 0.0 - 7th October 1996
C   VERSION 0.1 - 22nd October 1997 - major bug fix (+ some minor ones)
C   VERSION 0.1b - 3 August 1999, minor bug fixes & modifications by GPS
C
C - NEV IS THE NUMBER OF EVENTS TO GENERATE
C - S IS THE TOTAL LEPTON-HADRON CENTRE-OF-MASS ENERGY SQUARED
C - NFL IS THE NUMBER OF FLAVOURS
C - USER IS THE NAME OF ROUTINE TO ANALYSE THE EVENTS
C   IT TAKES SIX ARGUMENTS:
C     SUBROUTINE USER(N,NA,ITYPE,P,S,WEIGHT)
C     N = THE NUMBER OF PARTONS (FIRST IS INCOMING, N-1 OUTGOING)
C     NA = THE ORDER IN ALPHA
C     ITYPE = THE TYPE: 0=TREE LEVEL, 1=SUBTRACTION, 2=FINITE VIRTUAL,
C                       3=FINITE COLLINEAR
C     P(4,7) = THEIR MOMENTA (P(1-4,I) IS 4-MOMENTUM OF PARTICLE I,
C              P(1-4,5) IS TOTAL EXCHANGED MOMENTUM (IE FIRST-REST)
C              P(1-4,6/7) IS THE 4-MOM OF THE INCOMING/OUTGOING LEPTON)
C     S = THE TOTAL LEPTON-HADRON CENTRE-OF-MASS ENERGY SQUARED
C     WEIGHT(I) = THE WEIGHT IN THE INTEGRAL INITIATED BY PARTONS OF
C                 TYPE I (0=G,1=D,2=U,3=S,...,6=T,-1=DBAR,...,-6=TBAR)
C   IT IS ALSO CALLED AT THE END OF EACH FULL EVENT WITH N=0
C   A SIMPLE EXAMPLE (CALLED DEMO) IS GIVEN BELOW
C - CUTS IS THE NAME OF ROUTINE THAT PROVIDES CUTS ON X-Q**2-Y
C   IT TAKES ONE ARGUMENT AND PROVIDES SIX OUTPUTS:
C     SUBROUTINE CUTS(S,XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX)
C   IF XMIN=XMAX, THE RESULTING CROSS-SECTION IS DIFFERENTIAL IN X
C   IF XMAX>XMIN, IT IS INTEGRATED OVER THE GIVEN X RANGE
C   IF XMAX=0, THE FULL ALLOWED X RANGE IS INTEGRATED OVER
C   AND LIKEWISE FOR Q2 AND Y.
C   AT MOST TWO OF THE VARIABLES CAN BE FIXED, BUT ALL CAN HAVE LIMITS
C   A LOWER LIMIT ON Q2 MUST EITHER BE EXPLICITLY GIVEN,
C   OR IMPLIED BY THE OTHER TWO LIMITS
C
C   IT IS IMPORTANT TO THINK ABOUT THE EFFICIENCY OF THE USER ROUTINE,
C   AS IT IS CALLED 14 TIMES PER EVENT
C  (THIS IS 8 3-PARTON CALLS, 4 2-PARTON CALLS AND 2 1-PARTON CALLS)
C
C
C GPS ADDITIONS: ADD ARGS NPOW1,NPOW2,CUTOFF_IN,SCALE_in 
C                TO SET THINGS THAT ARE COMMONLY CHANGED
C
      INTEGER NEV,NFL,NPERM3,NPERM4,I,J
      PARAMETER (NPERM3=1,NPERM4=6)
      DOUBLE PRECISION S,CA_in, CF_in, TR_in,
     $     WTWO,WTHR,WFOR,JTHR,JFOR,JTMP,MTWO(-6:6),MTHR(-6:6),
     $     MFOR(-6:6),VTWO(-6:6),VTHR(-6:6),CTHR(-6:6),CFOR(-6:6),
     $     STHR(-6:6,NPERM3),SFOR(-6:6,NPERM4),WEIGHT(-6:6),ZERO(-6:6),
     $     P(4,7),Q(4,7,NPERM4),NRM,NPOW1,NPOW2,CUTOFF_IN,SCALE_in
      INTEGER SCHEME,NF,ORDER
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      CHARACTER*10 date, time, zone
      INTEGER(8) values(8)
!     AK To set the seeds by hand
      INTEGER USD1, USD2
      DOUBLE PRECISION USEEDS(2) 
      COMMON  /SAMPLE/ XPOW,NPOW
      EXTERNAL USER,CUTS
      DATA ZERO/13*0/
!     AK: Scale variations
      LOGICAL SCALE_VAR, SCALE_VAR_in
      DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
C---PRINT OPENING MESSAGE
      WRITE (6,'(/2A)')  ' This is DISENT, a program for calculating',
     $     ' jet quantities in'
      WRITE (6,'(2A)')   ' deep inelastic scattering to',
     $     ' next-to-leading order in alpha_s'
      WRITE (6,'(A)')    ' If you use this program, please reference:'
      WRITE (6,'(2A)')   ' S.Catani & M.H.Seymour,',
     $     ' Nucl. Phys. B485 (1997) 291'
      WRITE (6,'(/A)')  ' Written by Mike Seymour, August 1996'
      WRITE (6,'(A/)')  ' Version 0.1, October 1997'
      WRITE (6,'(A)') ' Minor modifications by Gavin Salam, August 1999'
      WRITE (6,'(A)')
     $     ' Minor modifications by Alexander Karlberg, March 2023'
      WRITE (6,'(A)') ' Including bug fix in SUBFOR as reported in'
      WRITE (6,'(A)') ' 2005.10705 and 2010.07354'
      WRITE (6,'(A,I10)') '    NEV=',NEV
C---INITIALIZE COLOUR FACTORS AND OTHER CONSTANTS
      CF=CF_in
      CA=CA_in
      TR=TR_in
      if (CA.ge.3.1.or.CA.lt.0d0.or.CF.ge.1.4d0.or.CF.lt.0d0) then
         write(0,*) 'Requested values of CF,CA,TR were'
         write(0,*) CF,CA,TR
         write(0,*) 'They do not make sense and will be reset!'
         CF = 4d0/3
         CA = 3d0
         TR = 0.5d0
      end if
      NF=NFL
      write(6,
     $   '(" CA =",f10.6,"  CF =",f10.6,"  TR =",f10.6,"  nf =",i3)')
     $     ca,cf,tr,nf
      PI=ATAN(1D0)*4
      PISQ=PI**2
      HF=0.5D0
      CUTOFF=CUTOFF_IN
      write(6,*) 'CUTOFF =',CUTOFF
      EQ(0)=0
      EQ(1)=-1D0/3
      EQ(2)=EQ(1)+1
      DO I=1,6
        IF (I.GT.2) EQ(I)=EQ(I-2)
        EQ(-I)=-EQ(I)
      ENDDO
C---SCHEME IS 0 FOR MSbar AND OTHERWISE FOR DIS
      SCHEME=0
      write(6,*) 'SCHEME =', SCHEME
C---SCALE IS FACTORIZATION SCALE**2/Q**2
c~       SCALE=2
      SCALE=SCALE_in
      SCALE_VAR = SCALE_VAR_in
!     AK: This converts to picobarn
      NRM=3.8937966d8
C---PARAMETERS RELATED TO IMPORTANCE SAMPLING
      NPOW(1)=NPOW1
      NPOW(2)=NPOW2
      write(6,*) 'NPOW =',NPOW(1),NPOW(2)
      XPOW(1)=1-1D0/NPOW(1)
      XPOW(2)=1-1D0/NPOW(2)
      write(6,*) 'ORDER =',ORDER
!     AK Set the seed
      USEEDS(1) = USD1
      USEEDS(2) = USD2
      CALL RANGEN(0,USEEDS)
!     AK: Below some modifications: - Removed some calls to
!     contributions that we do not need since they are already in the
!     structure functions. This leads to a ~10% reduction in runtime at
!     NNLO and ~30% at NLO.
C---  START MAIN LOOP
      DO I=1,NEV
         IF (MOD(I,100 000).EQ.1) CALL RANGEN(-I,SFOR)
C---  GENERATE A TWO-PARTON STATE
         SCL_WEIGHT = 1D0 
         CALL GENTWO(S,P,WTWO,CUTS,*1000)
C---  EVALUATE THE TWO-PARTON TREE-LEVEL MATRIX ELEMENT
!         CALL MATTWO(P,MTWO)
C---  GIVE IT TO THE USER
!         CALL VECMUL(13,NRM*WTWO/NEV,MTWO,WEIGHT)
!         CALL USER(2,0,0,P,S,WEIGHT)
         
         IF(ORDER.GE.1) THEN
C---  EVALUATE THE TWO-PARTON ONE-LOOP MATRIX ELEMENT
            CALL VIRTWO(S,P,VTWO,*1000)
C---  GIVE IT TO THE USER
!            CALL VECMUL(13,NRM*WTWO/NEV,VTWO,WEIGHT)
!            CALL USER(2,1,2,P,S,WEIGHT)
C---  EVALUATE THE THREE-PARTON COLLINEAR SUBTRACTION
            CALL COLTHR(S,P,CTHR,*1000)
C---  GIVE IT TO THE USER
!            CALL VECMUL(13,NRM*WTWO/NEV,CTHR,WEIGHT)
!            CALL USER(3,1,3,P,S,WEIGHT)
C---  GENERATE A THREE-PARTON STATE
            CALL GENTHR(P,WTHR,*1000)
C---  CALCULATE THE JACOBIAN FACTOR AND SUBTRACTION CONFIGURATIONS
            JTHR=0
            DO J=1,NPERM3
               CALL SUBTHR(J,S,P,Q(1,1,J),STHR(-6,J),JTMP,*1000)
               JTHR=JTHR+JTMP
            ENDDO
C---  EVALUATE THE THREE-PARTON TREE-LEVEL MATRIX ELEMENT
            CALL MATTHR(P,MTHR)
C---  GIVE IT TO THE USER
            CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,MTHR,WEIGHT)
            CALL USER(3,1,0,P,S,WEIGHT)
C---  GIVE THE SUBTRACTION CONFIGURATIONS TO THE USER
!            DO J=1,NPERM3
!               CALL VECMUL(13,-NRM*WTWO*WTHR/JTHR/NEV,STHR(-6,J),WEIGHT)
!               CALL USER(3,1,1,Q(1,1,J),S,WEIGHT)
!            ENDDO
         ENDIF
C---  GPS MODIFICATION --------------------
         IF(ORDER.GE.2) THEN
C---  EVALUATE THE THREE-PARTON ONE-LOOP MATRIX ELEMENT
            CALL VIRTHR(S,P,VTHR,*1000)
C---  GIVE IT TO THE USER
            CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,VTHR,WEIGHT)
            CALL USER(3,2,2,P,S,WEIGHT)
C---  EVALUATE THE FOUR-PARTON COLLINEAR SUBTRACTION
            CALL COLFOR(S,P,CFOR,*1000)
C---  GIVE IT TO THE USER
            CALL VECMUL(13,NRM*WTWO*WTHR/JTHR/NEV,CFOR,WEIGHT)
            CALL USER(4,2,3,P,S,WEIGHT)
C---  GENERATE A FOUR-PARTON STATE
            SCL_WEIGHT = 1D0 
            CALL GENFOR(P,WFOR,*1000)
C---  CALCULATE THE JACOBIAN FACTOR AND SUBTRACTION CONFIGURATIONS
            JFOR=0
            DO J=1,NPERM4
               CALL SUBFOR(J,S,P,Q(1,1,J),SFOR(-6,J),JTMP,*1000)
               JFOR=JFOR+JTMP
            ENDDO
C---  EVALUATE THE FOUR-PARTON TREE-LEVEL MATRIX ELEMENT
            CALL MATFOR(P,MFOR)
C---  GIVE IT TO THE USER
            CALL VECMUL(13,NRM*WTWO*WTHR*WFOR/JFOR/NEV,MFOR,WEIGHT)
            CALL USER(4,2,0,P,S,WEIGHT)
C---  GIVE THE SUBTRACTION CONFIGURATIONS TO THE USER
            DO J=1,NPERM4
               CALL VECMUL(13,
     $              -NRM*WTWO*WTHR*WFOR/JFOR/NEV,SFOR(-6,J),WEIGHT)
               CALL USER(4,2,1,Q(1,1,J),S,WEIGHT)
            ENDDO
         END IF
C---  TELL THE USER THAT THE EVENT IS COMPLETE
 1000    CALL USER(0,0,0,P,S,ZERO)
      ENDDO
      END
C-----------------------------------------------------------------------
      FUNCTION DOT(P,I,J)
      IMPLICIT NONE
C---RETURN THE DOT PRODUCT OF P(*,I) AND P(*,J)
      INTEGER I,J
      DOUBLE PRECISION DOT,P(4,7)
      DOT=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
      END
C-----------------------------------------------------------------------
      SUBROUTINE VECMUL(N,S,A,B)
      IMPLICIT NONE
C---MULTIPLY THE VECTOR A BY THE SCALAR S TO GIVE THE VECTOR B
      INTEGER N,I
      DOUBLE PRECISION S,A(N),B(N)
      DO I=1,N
        B(I)=S*A(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENTWO(S,P,W,CUTS,*)
      IMPLICIT NONE
C---GENERATE A TWO-PARTON CONFIGURATION
C   WITHIN THE PHASE-SPACE LIMITS GIVEN IN CUTS()
C   THE WEIGHT GIVES THE TOTAL VOLUME OF PHASE-SPACE
      DOUBLE PRECISION S,P(4,7),W,R(2),XMIN,XMAX,QMIN,QMAX,YMIN,YMAX,
     $     Q2,QJAC,Y,YJAC,E
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      EXTERNAL CUTS
      CALL RANGEN(2,R)
C---FIRST GET THE USER'S CUTS
      CALL CUTS(S,XMIN,XMAX,QMIN,QMAX,YMIN,YMAX)
C---CHECK THAT SYSTEM IS NOT 0VER-CONSTRAINED
      IF (XMIN.EQ.XMAX.AND.XMIN.NE.0.AND.
     $    YMIN.EQ.YMAX.AND.YMIN.NE.0.AND.
     $    QMIN.EQ.QMAX.AND.QMIN.NE.0)
     $     STOP 'AT MOST TWO KINEMATIC VARIABLES CAN BE CONSTRAINED!'
C---INSERT KINEMATIC LIMITS WHERE NONE ARE GIVEN
      IF (XMAX.LE.0) XMAX=1
      IF (YMAX.LE.0) YMAX=1
      IF (QMAX.LE.0) QMAX=S
C---CHECK IF Q2 IS CONSTRAINED BY X AND Y
      IF (XMAX*YMAX*S.LT.QMAX) QMAX=XMAX*YMAX*S
      IF (XMIN*YMIN*S.GT.QMIN) QMIN=XMIN*YMIN*S
C---CHECK THERE IS SOME RANGE AVAILABLE IN Q2
      IF (QMAX.LT.QMIN) STOP 'NO PHSP AVAILABLE (QMAX.LT.QMIN)'
      IF (QMIN.EQ.0) STOP 'Q2 MUST BE CONSTRAINED LARGER THAN ZERO!'
C---GENERATE A Q2 VALUE
      IF (QMAX.EQ.QMIN) THEN
        Q2=QMAX
        QJAC=1
      ELSE
        Q2=(QMAX/QMIN)**R(1)*QMIN
        QJAC=LOG(QMAX/QMIN)*Q2
      ENDIF
C---CHECK IF Y IS CONSTRAINED BY X AND Q2
      IF (XMIN*YMAX*S.GT.Q2) YMAX=Q2/(XMIN*S)
      IF (XMAX*YMIN*S.LT.Q2) YMIN=Q2/(XMAX*S)
C---CHECK THERE IS SOME RANGE AVAILABLE IN Y
      IF (YMAX.LT.YMIN) STOP 'NO PHSP AVAILABLE (YMAX.LT.YMIN)'
C---GENERATE A Y VALUE
      IF (YMAX.EQ.YMIN) THEN
        Y=YMAX
        YJAC=1
      ELSE
        Y=(YMAX/YMIN)**R(2)*YMIN
        YJAC=LOG(YMAX/YMIN)*Y
      ENDIF
C---  CONSTRUCT MOMENTA (IN THE BREIT FRAME, AS IT HAPPENS)
      P = 0
      E=SQRT(Q2)/2
      P(1,1)=   0
      P(2,1)=   0
      P(3,1)=   E
      P(4,1)=   E
      P(1,2)=   0
      P(2,2)=   0
      P(3,2)=  -E
      P(4,2)=   E
      P(1,5)=   0
      P(2,5)=   0
      P(3,5)=-2*E
      P(4,5)=   0
      P(1,6)=   E/Y*2*SQRT(1-Y)
      P(2,6)=   0
      P(3,6)=  -E
      P(4,6)=   E/Y*(2-Y)
      P(1,7)=   E/Y*2*SQRT(1-Y)
      P(2,7)=   0
      P(3,7)=   E
      P(4,7)=   E/Y*(2-Y)

C---CALCULATE WEIGHT
      W=Y/(16*PI*Q2**2)*QJAC*YJAC
C---  AB-return x*dsigma/dx, same as NLOJET++
      if (xmax==xmin) w=w*y
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENTHR(P,W,*)
      IMPLICIT NONE
C---GENERATE A THREE-PARTON CONFIGURATION FROM A GIVEN TWO-PARTON ONE
      DOUBLE PRECISION P(4,7),W,R(2),X,XJAC,Z,EMSQ
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      CALL RANGEN(2,R)
C---FIND OUT WHAT X VALUE WAS GENERATED AND THE JACOBIAN FACTOR
      CALL GETCOL(X,XJAC)
C---GENERATE A Z VALUE
      Z=1-R(1)**NPOW(1)
      IF (R(2).GT.0.5) Z=1-Z
C---  ENFORCE INVARIANT MASS CUTOFF
      IF (Z.LT.CUTOFF.OR.1-Z.LT.CUTOFF) RETURN 1
C---  GENERATE THEIR MOMENTA
      CALL GENDEC(P,2,3,Z,*999)
C---CALCULATE WEIGHT
      EMSQ=P(4,5)**2-P(3,5)**2-P(2,5)**2-P(1,5)**2
      W=-EMSQ/(16*PISQ)
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENFOR(P,W,*)
      IMPLICIT NONE
C---GENERATE A FOUR-PARTON CONFIGURATION FROM A GIVEN THREE-PARTON ONE
      INTEGER EMIT,OMIT,I
      DOUBLE PRECISION P(4,7),W,R(6),X,XJAC,Z,EMSQ,XX(2),T
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      CALL RANGEN(6,R)
C---CHOOSE WHICH PARTON TO SPLIT UNIFORMLY
      EMIT=INT(R(1)*2)+2
C---AND WHICH OF THE OTHERS TO BALANCE ITS MOMENTUM WITH
      OMIT=INT(R(2)*2)+1
      IF (OMIT.EQ.EMIT) OMIT=3
C---IF INCOMING PARTON IS THE SPECTATOR GENERATE AN INITIAL-STATE DIPOLE
      IF (OMIT.EQ.1) THEN
C---FIND OUT WHAT X VALUE WAS GENERATED AND THE JACOBIAN FACTOR
        CALL GETCOL(X,XJAC)
C---GENERATE A Z VALUE
        Z=1-R(3)**NPOW(2)
        IF (R(4).GT.0.5) Z=1-Z
C---ENFORCE INVARIANT MASS CUTOFF
        IF (Z.LT.CUTOFF.OR.1-Z.LT.CUTOFF) RETURN 1
C---GENERATE THEIR MOMENTA
        CALL GENDEC(P,EMIT,4,Z,*999)
C---OTHERWISE GENERATE A FINAL-STATE DIPOLE
      ELSE
C---AND FORGET THE PROPOSED COLLINEAR SPLITTING
        DO I=1,4
          P(I,1)=P(I,1)-P(I,4)
        ENDDO
C---GENERATE X1,X2 VALUES
        XX(1)=1-MIN(R(3),R(4))**NPOW(2)
        XX(2)=2-XX(1)-(R(5)+(1-R(5))*(1-XX(1))**(1D0/NPOW(2)))**NPOW(2)
C---ENFORCE INVARIANT MASS CUTOFFS
        IF (1-XX(1).LT.CUTOFF.OR.1-XX(2).LT.CUTOFF
     $       .OR.XX(1)+XX(2)-1.LT.CUTOFF) RETURN 1
C---GENERATE THEIR MOMENTA
        CALL GENDIP(P,EMIT,OMIT,4,XX,*999)
      ENDIF
C---HALF THE TIME, SWAP PARTICLES 3 AND 4 SO THEY ARE SYMMETRIC
      IF (R(6).GT.0.5) THEN
        DO I=1,4
          T=P(I,3)
          P(I,3)=P(I,4)
          P(I,4)=T
        ENDDO
      ENDIF
C---CALCULATE WEIGHT
      EMSQ=P(4,5)**2-P(3,5)**2-P(2,5)**2-P(1,5)**2
      W=-EMSQ/(16*PISQ)
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENDEC(P,I1,I2,Z,*)
      IMPLICIT NONE
C---GENERATE A 2->2 DECAY
      DOUBLE PRECISION P(4,7),Z,R(1),Q(4,2),PTDQ,C(4),D(4)
      INTEGER I1,I2,I
      LOGICAL FIRST
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DATA FIRST/.TRUE./
      CALL RANGEN(1,R)
      PTDQ=Z*(1-Z)
      IF (PTDQ.LT.0) THEN
        IF (FIRST) THEN
          FIRST=.FALSE.
          WRITE (*,*) 'Numerical error in GENDEC!'
          WRITE (*,*) 'Please report this to seymour@surya11.cern.ch,'
          WRITE (*,*) 'giving the date of this version,'
          WRITE (*,*) 'and the latest value of ISEED.'
          WRITE (*,*) 'Event generation has not been affected.'
        ELSE
          WRITE (*,*) 'Another numerical error in GENDEC!'
        ENDIF
        RETURN 1
      ENDIF
      PTDQ=SQRT(PTDQ)
      CALL GTPERP(PTDQ,P,I2,I1,6,C,D)
      DO I=1,4
        Q(I,1)=Z*P(I,I1)+(1-Z)*P(I,I2)
     $       +COS(R(1)*2*PI)*C(I)+SIN(R(1)*2*PI)*D(I)
        Q(I,2)=(1-Z)*P(I,I1)+Z*P(I,I2)
     $       -COS(R(1)*2*PI)*C(I)-SIN(R(1)*2*PI)*D(I)
      ENDDO
      DO I=1,4
        P(I,I1)=Q(I,1)
        P(I,I2)=Q(I,2)
      ENDDO
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENDIP(P,EMIT,OMIT,THRD,X,*)
      IMPLICIT NONE
C---GENERATE A 2->3 DIPOLE EMISSION
      DOUBLE PRECISION P(4,7),X(2),R(1),Q(4,2),PTDQ,C(4),D(4)
      INTEGER EMIT,OMIT,THRD,I
      LOGICAL FIRST
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DATA FIRST/.TRUE./
      CALL RANGEN(1,R)
      PTDQ=(X(1)+X(2)-1)*(1-X(1))*(1-X(2))/X(1)**2
      IF (PTDQ.LT.0) THEN
        IF (FIRST) THEN
          FIRST=.FALSE.
          WRITE (*,*) 'Numerical error in GENDIP!'
          WRITE (*,*) 'Please report this to seymour@surya11.cern.ch,'
          WRITE (*,*) 'giving the date of this version,'
          WRITE (*,*) 'and the latest value of ISEED.'
          WRITE (*,*) 'Event generation has not been affected.'
        ELSE
          WRITE (*,*) 'Another numerical error in GENDIP!'
        ENDIF
        RETURN 1
      ENDIF
      PTDQ=SQRT(PTDQ)
      CALL GTPERP(PTDQ,P,EMIT,OMIT,6,C,D)
      DO I=1,4
        Q(I,1)=X(1)*P(I,OMIT)
        Q(I,2)=(1-(1-X(2))/X(1))*P(I,EMIT)
     $       +(1-X(1))*(1-X(2))/X(1)*P(I,OMIT)
     $       +COS(R(1)*2*PI)*C(I)+SIN(R(1)*2*PI)*D(I)
      ENDDO
      DO I=1,4
        P(I,THRD)=P(I,EMIT)+P(I,OMIT)-Q(I,1)-Q(I,2)
        P(I,OMIT)=Q(I,1)
        P(I,EMIT)=Q(I,2)
      ENDDO
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE GENCOL(I,X,XJAC,XMIN)
      IMPLICIT NONE
C---GENERATE AN X VALUE AND STORE IT FOR LATER RETRIEVAL
      INTEGER I
      DOUBLE PRECISION X,XJAC,XMIN,XL,XLJAC,R(2)
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      SAVE XL,XLJAC
      CALL RANGEN(2,R)
      IF (R(1).LT.0.5) THEN
        X=1-(1-XMIN)*R(2)**NPOW(I)
      ELSE
        X=XMIN**R(2)
      ENDIF
      IF (X.EQ.1) THEN
        XJAC=0
      ELSE
        XJAC=1/(0.5/(-X*LOG(XMIN))
     $       +0.5*((1-XMIN)/(1-X))**XPOW(I)/(NPOW(I)*(1-XMIN)))
      ENDIF
      XL=X
      XLJAC=XJAC
      RETURN
      ENTRY GETCOL(X,XJAC)
C---RETURN THE SAME VALUES AS LAST TIME
      X=XL
      XJAC=XLJAC
      END
C-----------------------------------------------------------------------
      SUBROUTINE GTPERP(PTDQ,P,I,J,K,C,D)
      IMPLICIT NONE
C---FIND THE VECTORS PERPENDICULAR TO P(I) AND P(J)
C   C AND D ARE PURELY SPACE-LIKE VECTORS IN THE P(I)+P(J) CMF,
C   WITH C IN THE SAME PLANE AS P(K) AND D PERPENDICULAR TO IT,
C   BOTH HAVING LENGTH PTDQ*SQRT(2*DOT(P,I,J))
      DOUBLE PRECISION PTDQ,P(4,7),C(4),D(4),PTF,DIJ,DIK,DJK,DOT,EPS4
      INTEGER I,J,K,L
      DIJ=DOT(P,I,J)
      DIK=DOT(P,I,K)
      DJK=DOT(P,J,K)
      PTF=PTDQ/SQRT(DIK*DJK)
      DO L=1,4
        C(L)=PTF*(DIJ*P(L,K)-DJK*P(L,I)-DIK*P(L,J))
      ENDDO
      DO L=1,4
        D(L)=EPS4(L,P(1,I),P(1,J),C)/DIJ
      ENDDO
      END
C-----------------------------------------------------------------------
      FUNCTION EPS4(I,A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION EPS4,EPS3,A(4),B(4),C(4),AA(3),BB(3),CC(3)
      INTEGER I,J,K,S(4)
      DATA S/+1,-1,+1,+1/
      J=1
      DO K=1,3
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        CC(K)=C(J)
        J=J+1
      ENDDO
      EPS4=0
      DO J=1,3
        EPS4=EPS4+CC(J)*EPS3(J,AA,BB)
      ENDDO
      EPS4=S(I)*EPS4
      END
C-----------------------------------------------------------------------
      FUNCTION EPS3(I,A,B)
      IMPLICIT NONE
      DOUBLE PRECISION EPS3,A(3),B(3),AA(2),BB(2)
      INTEGER I,J,K,S(3)
      DATA S/+1,-1,+1/
      J=1
      DO K=1,2
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        J=J+1
      ENDDO
      EPS3=S(I)*(AA(1)*BB(2)-AA(2)*BB(1))
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATTWO(P,M)
      IMPLICIT NONE
C---EVALUATE THE TWO-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I
      DOUBLE PRECISION P(4,7),M(-6:6),Q,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      Q=4*(4*PI/137)**2/DOT(P,5,5)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
      DO I=-6,6
        M(I)=EQ(I)**2*Q
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATTHR(P,M)
      IMPLICIT NONE
C---EVALUATE THE THREE-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I
      DOUBLE PRECISION P(4,7),M(-6:6),QQ,GQ,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      QQ=8*(4*PI/137)**2*
     $     (DOT(P,1,6)**2+DOT(P,1,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*CF/(-4*DOT(P,2,3)*DOT(P,1,3)*DOT(P,5,5))
      GQ=8*(4*PI/137)**2*
     $     (DOT(P,3,6)**2+DOT(P,3,7)**2+DOT(P,2,7)**2+DOT(P,2,6)**2)
     $     *16*PISQ*TR/(-4*DOT(P,2,1)*DOT(P,3,1)*DOT(P,5,5))
      DO I=-6,6
        M(I)=EQ(I)**2*QQ
      ENDDO
      M(0)=0
      DO I=1,NF
        M(0)=M(0)+EQ(I)**2*GQ
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE CONTHR(P,V,VV,Q,QBAR,G,M)
      IMPLICIT NONE
C---EVALUATE THE CONTRACTION OF THE THREE-PARTON MATRIX ELEMENT SQUARED
C   FOR THE GIVEN CONFIGURATION WITH THE VECTOR V.  IT IS ASSUMED THAT
C   V IS PERPENDICULAR TO THE GLUON, V.P3.EQ.0, AND THAT V.V.NE.0 .
C   THE NORMALIZATION IS SUCH THAT THE AVERAGE OVER V'S AZIMUTH AROUND
C   P3 IS EQUAL TO HALF OF MATTHR
C
C GPS: added VV as one of the arguments, to allow a better calculation 
C      of it by the calling routines
      INTEGER I,J,Q,QBAR,G
      DOUBLE PRECISION P(4,7),V(4),M,DOT,VDOT,EMSQ,OMSQ,VV,Y1,Y2,L1,L2,
     $     DOTS,T1,T2,T3,T4,T,V1,V2,V6
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      VDOT(I)=P(4,I)*V(4)-P(3,I)*V(3)-P(2,I)*V(2)-P(1,I)*V(1)
      EMSQ=DOT(P,5,5)
      OMSQ=1/EMSQ
C---  let VV be calculated by the calling routine
C      VV=V(4)**2-V(3)**2-V(2)**2-V(1)**2
      Y1=2*DOT(P,ABS(QBAR),ABS(G))*SIGN(1,QBAR*G)
      Y2=2*DOT(P,ABS(Q),ABS(G))*SIGN(1,Q*G)
      L1=2*DOT(P,ABS(Q),6)*SIGN(1,Q)
      L2=2*DOT(P,ABS(QBAR),6)*SIGN(1,QBAR)
      V1=VDOT(ABS(Q))*SIGN(1,Q)
      V2=VDOT(ABS(QBAR))*SIGN(1,QBAR)
      V6=VDOT(6)
      DOTS=1/(VV*(Y1+Y2)**2)
      T=CF*16*PISQ/(Y1*Y2)**2
      T1=T*4*(Y1**2*(L2**2+(EMSQ-L2)**2)+Y2**2*(L1**2+(EMSQ-L1)**2)
     $     +2*Y1*Y2*(L1*EMSQ+L2*EMSQ-2*L1*L2))*OMSQ**2*DOTS
      T2=T*4*(Y1*(EMSQ-2*L2)-Y2*(EMSQ-2*L1))*Y1*Y2*OMSQ**2*DOTS
      T3=T*8*(Y1*Y2*OMSQ)**2*DOTS
      T4=T*((EMSQ-Y1-L1-L2)**2+(EMSQ-Y2-L1-L2)**2)
     $     *Y1*Y2*OMSQ
      M=-T1*(V1*Y1-V2*Y2)**2
     $  -T2*2*(V1*Y1-V2*Y2)*(-V6*(Y2+Y1)+(EMSQ-L1-L2)*(V1+V2))
     $  -T3*(-V6*(Y2+Y1)+(EMSQ-L1-L2)*(V1+V2))**2
     $  +T4
      m=m*(4*pi/137)**2
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATFOR(P,M)
      IMPLICIT NONE
C---EVALUATE THE FOUR-PARTON MATRIX ELEMENT SQUARED FOR THE GIVEN
C   CONFIGURATION.
      INTEGER I,J
      DOUBLE PRECISION P(4,7),M(-6:6),ERTA,ERTB,ERTC,ERTD,ERTE,
     $     LEIA,LEIB,LEIC,LEID,LEIE,
     $     A,B,C,DS,D1,D2,E,Q,G,QQ,EMSQ,DOT
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      EMSQ=-DOT(P,5,5)
      A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4)
     $    +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)
     $    +ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
      B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4)
     $    +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)
     $    +ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))
      C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4)
     $    +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*(
     $     ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)
     $    +ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
      Q=HF*(CF*A+(CF-CA/2)*B+CA*C)
      A=2*(LEIA(P,P(1,6),2,3,-1,4)+LEIA(P,P(1,6),3,2,-1,4)
     $    +LEIA(P,P(1,6),2,3,4,-1)+LEIA(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTA(P,2,3,-1,4)+ERTA(P,3,2,-1,4)
     $    +ERTA(P,2,3,4,-1)+ERTA(P,3,2,4,-1))
      B=2*(LEIB(P,P(1,6),2,3,-1,4)+LEIB(P,P(1,6),3,2,-1,4)
     $    +LEIB(P,P(1,6),2,3,4,-1)+LEIB(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTB(P,2,3,-1,4)+ERTB(P,3,2,-1,4)
     $    +ERTB(P,2,3,4,-1)+ERTB(P,3,2,4,-1))
      C=2*(LEIC(P,P(1,6),2,3,-1,4)+LEIC(P,P(1,6),3,2,-1,4)
     $    +LEIC(P,P(1,6),2,3,4,-1)+LEIC(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $     ERTC(P,2,3,-1,4)+ERTC(P,3,2,-1,4)
     $    +ERTC(P,2,3,4,-1)+ERTC(P,3,2,4,-1))
      G=-(CF*A+(CF-CA/2)*B+CA*C)
      D1=2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*(
     $      ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))
      Q=Q+NF*TR*D1
      D2=2*(LEID(P,P(1,6),-1,4,2,3)+LEID(P,P(1,6),2,3,-1,4))-EMSQ/2*(
     $      ERTD(P,-1,4,2,3)+ERTD(P,2,3,-1,4))
      QQ=TR*D2
c$$$      DS=2*(LEID(P,P(1,6),4,-1,2,3)+LEID(P,P(1,6),2,3,4,-1)
c$$$     $     +LEID(P,P(1,6),-1,4,3,2)+LEID(P,P(1,6),3,2,-1,4))-EMSQ/2*(
c$$$     $      ERTD(P,4,-1,2,3)+ERTD(P,2,3,4,-1)
c$$$     $     +ERTD(P,-1,4,3,2)+ERTD(P,3,2,-1,4))
c$$$      Q=Q+HF*TR*(DS-D1-D2)
      E=2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1)
     $    +LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1)
     $    +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4)
     $    +LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4))-EMSQ/2*(
     $     ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)
     $    +ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1)
     $    +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)
     $    +ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
      Q=Q+HF*(CF-CA/2)*E
C---INCLUDE EXTERNAL FACTORS
      Q=Q*256*PI**4*CF/EMSQ
      Q=Q*(4*PI/137)**2*4/EMSQ
      G=G*256*PI**4*TR/EMSQ
      G=G*(4*PI/137)**2*4/EMSQ
      QQ=QQ*256*PI**4*CF/EMSQ
      QQ=QQ*(4*PI/137)**2*4/EMSQ
      DO I=-6,6
        M(I)=EQ(I)**2*Q
        DO J=1,NF
          M(I)=M(I)+EQ(J)**2*QQ
        ENDDO
      ENDDO
      M(0)=0
      DO I=1,NF
        M(0)=M(0)+EQ(I)**2*G
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VIRTWO(S,P,V,*)
      IMPLICIT NONE
C---CALCULATE THE TWO-PARTON MATRIX-ELEMENT AT NEXT-TO-LEADING ORDER
      INTEGER I
      DOUBLE PRECISION S,P(4,7),V(-6:6),M(-6:6),O,X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,DOT
      PARAMETER (O=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
            DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      DOUBLE PRECISION QQscl(3) ,GQscl(3),QGscl(3) ,GGscl(3)
      LOGICAL SCALE_VAR
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTWO(P,M)
C---SUM OF FACTORIZING VIRTUAL CROSS-SECTION AND SUBTRACTION COUNTERTERM
      QQ=CF*(2-PISQ)
      GG=0
C---THE NON-FACTORIZING VIRTUAL CROSS-SECTION
      DO I=-6,6
        V(I)=0
      ENDDO
C---GENERATE A COLLINEAR EMISSION
      XMIN=2*DOT(P,1,6)/S
      CALL GENCOL(1,X,XJAC,XMIN)
C---ENFORCE INVARIANT MASS CUTOFF
      IF (1-X.LT.CUTOFF) RETURN 1
      RETURN ! AK: We just need the phase space for the next routines
C---CALCULATE THE COLLINEAR COUNTERTERM
      GQ=0
      QG=0
      KQF=1.5
      if(SCALE_VAR) then
         QQscl = QQ
         GQscl = GQ
         QGscl = QG
         GGscl = GG
         do i = 1,3
            SCL_WEIGHT(i,:) = V(:)
         enddo
         CALL KPFUNS_SCL_VAR(-X,XJAC,XMIN,KQF,O,O,O,QQscl,GQscl,QGscl
     $        ,GGscl)
C---  THE TOTAL
         SCL_WEIGHT(:,0)=SCL_WEIGHT(:,0)+GGscl(:)*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               SCL_WEIGHT(:,I)=SCL_WEIGHT(:,I)+QGscl(:)*M(0)+QQscl(:)
     $              *M(I)
               IF (ABS(I).LE.NF) SCL_WEIGHT(:,0)=SCL_WEIGHT(:,0)
     $              +GQscl(:)*M(I)
            ENDIF
         ENDDO
         V(:) = SCL_WEIGHT(1,:) 
         do i = 1,3
            SCL_WEIGHT(i,:) = SCL_WEIGHT(i,:) / V(:) 
         enddo
      else
         CALL KPFUNS(-X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---  THE TOTAL
         V(0)=V(0)+GG*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               V(I)=V(I)+QG*M(0)+QQ*M(I)
               IF (ABS(I).LE.NF) V(0)=V(0)+GQ*M(I)
            ENDIF
         ENDDO
      endif
!      CALL KPFUNS(-X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---THE TOTAL
!      V(0)=V(0)+GG*M(0)
!      DO I=-6,6
!        IF (I.NE.0) THEN
!          V(I)=V(I)+QG*M(0)+QQ*M(I)
!          IF (ABS(I).LE.NF) V(0)=V(0)+GQ*M(I)
!        ENDIF
!      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VIRTHR(S,P,V,*)
      IMPLICIT NONE
C---CALCULATE THE THREE-PARTON MATRIX-ELEMENT AT NEXT-TO-LEADING ORDER
      INTEGER I
      DOUBLE PRECISION S,P(4,7),V(-6:6),M(-6:6),X,XJAC,XMIN, QQ,GQ,QG,GG
     $     ,KQF,KGF,PQF,PGF,L12,L13,L23,DOT,ERTV,LEIV,EMSQ
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      DOUBLE PRECISION QQscl(3) ,GQscl(3),QGscl(3) ,GGscl(3)
      LOGICAL SCALE_VAR
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
      EMSQ=-DOT(P,5,5)
      L12=LOG(2*DOT(P,1,2)/EMSQ)
      L13=LOG(2*DOT(P,1,3)/EMSQ)
      L23=LOG(2*DOT(P,2,3)/EMSQ)
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTHR(P,M)
C---THE NON-FACTORIZING VIRTUAL CROSS-SECTION
      QQ=-((4*PI/137)**2*4/EMSQ)*
     $     (2*LEIV(P,P(1,6),2,-1,3)-EMSQ/2*ERTV(P,2,-1,3))
      GG=TR/CF*((4*PI/137)**2*4/EMSQ)*
     $     (2*LEIV(P,P(1,6),2,3,-1)-EMSQ/2*ERTV(P,2,3,-1))
      DO I=-6,6
        V(I)=EQ(I)**2*QQ
      ENDDO
      V(0)=0
      DO I=1,NF
        V(0)=V(0)+EQ(I)**2*GG
      ENDDO
C---SUM OF FACTORIZING VIRTUAL CROSS-SECTION AND SUBTRACTION COUNTERTERM
      QQ=CF*2+CA*50D0/9-TR*NF*16D0/9-CF*PISQ
     $     -3*(CF-CA/2)*L12-(5*CA-TR*NF)/3*(L13+L23)
      GG=CF*2+CA*50D0/9-TR*NF*16D0/9-CA*PISQ
     $     -3*(CF-CA/2)*L23-(5*CA-TR*NF)/3*(L12+L13)
C---GENERATE A COLLINEAR EMISSION
      XMIN=2*DOT(P,1,6)/S
      CALL GENCOL(2,X,XJAC,XMIN)
C---ENFORCE INVARIANT MASS CUTOFF
      IF (1-X.LT.CUTOFF) RETURN 1
C---CALCULATE THE COLLINEAR COUNTERTERM
      GQ=0
      QG=0
      KQF=(1.5*(CF-CA/2)+0.5*(11D0/6*CA-2D0/3*NF*TR))/CF
      KGF=1.5
      PQF=-((CF-CA/2)*L12+CA/2*L13)/CF
      PGF=-(L12+L13)/2
      if(SCALE_VAR) then
         QQscl = QQ
         GQscl = GQ
         QGscl = QG
         GGscl = GG
         do i = 1,3
            SCL_WEIGHT(i,:) = V(:)
         enddo
         CALL KPFUNS_SCL_VAR(-X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQscl,GQscl
     $        ,QGscl,GGscl)
C---  THE TOTAL
         SCL_WEIGHT(:,0)=SCL_WEIGHT(:,0)+GGscl(:)*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               SCL_WEIGHT(:,I)=SCL_WEIGHT(:,I)+QGscl(:)*M(0)+QQscl(:)
     $              *M(I)
               IF (ABS(I).LE.NF) SCL_WEIGHT(:,0)=SCL_WEIGHT(:,0)
     $              +GQscl(:)*M(I)
            ENDIF
         ENDDO
         V(:) = SCL_WEIGHT(1,:) 
         do i = 1,3
            SCL_WEIGHT(i,:) = SCL_WEIGHT(i,:) / V(:) 
         enddo
      else
         CALL KPFUNS(-X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
C---  THE TOTAL
         V(0)=V(0)+GG*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               V(I)=V(I)+QG*M(0)+QQ*M(I)
               IF (ABS(I).LE.NF) V(0)=V(0)+GQ*M(I)
            ENDIF
         ENDDO
      endif
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ERTV(P,I,J,K)
      IMPLICIT NONE
C---RETURN THE ERT F FUNCTION FOR VIRTUAL TERMS THAT ARE NOT
C   TRIVIALLY PROPORTIONAL TO TREE-LEVEL
      INTEGER I,J,K
      DOUBLE PRECISION P(4,7),DOT,R,RR,X,Y,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,R1312,R2312,R2313
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))
      R1312=R(Y13,Y12)
      R2312=R(Y23,Y12)
      R2313=R(Y23,Y13)
      CFSUB=Y12/(Y12+Y13)+Y12/(Y12+Y23)+(Y12+Y23)/Y13+(Y12+Y13)/Y23
     $     +L13*(4*Y12**2+2*Y12*Y13+4*Y12*Y23+Y13*Y23)/(Y12+Y23)**2
     $     +L23*(4*Y12**2+2*Y12*Y23+4*Y12*Y13+Y13*Y23)/(Y12+Y13)**2
     $     -2*((Y12**2+(Y12+Y13)**2)/(Y13*Y23)*R2312
     $        +(Y12**2+(Y12+Y23)**2)/(Y13*Y23)*R1312
     $        +(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))
     $     -2*L12*(Y12**2/(Y13+Y23)**2+2*Y12/(Y13+Y23)))
      CASUB=L13*Y13/(Y12+Y23)+L23*Y23/(Y12+Y13)
     $     +((Y12**2+(Y12+Y13)**2)/(Y13*Y23)*R2312
     $     +(Y12**2+(Y12+Y23)**2)/(Y13*Y23)*R1312
     $     +(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))
     $     -2*L12*(Y12**2/(Y13+Y23)**2+2*Y12/(Y13+Y23)))
     $     -R2313*(Y13/Y23+Y23/Y13+2*Y12/(Y13*Y23))
      TFSUB=0
      ERTV=CF*(CF*CFSUB+CA*CASUB)
      CFSUB=CFSUB/ERTV
      CASUB=CASUB/ERTV
      ERTV=16*PISQ/EMSQ*ERTV
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LEIV(P,Q,I,J,K)
      IMPLICIT NONE
C---RETURN THE PART OF THE ONE-LOOP HADRONIC TENSOR NOT TRIVIALLY
C   PROPORTIONAL TO TREE-LEVEL ONE CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
      INTEGER I,J,K
      DOUBLE PRECISION P(4,7),Q(4),DOT,R,RR,X,Y,Z,DILOG,EMSQ,
     $     Y12,Y13,Y23,L12,L13,L23,A,B,C,D,P1Q,P2Q,P3Q,QQ,
     $     CF1,CF2,CF3,CF4,CA1,CA2,CA3,CA4,R1312,R2312,R2313
      PARAMETER (Z=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      COMMON  /SUBCOM/ CFSUB,CASUB,TFSUB,GGSUB,QQSUB,QPSUB
      RR(X,Y)=DILOG(X)+DILOG(Y)-DILOG(X*Y)-PISQ/6
      R(X,Y)=RR((1-X)/Y,(1-Y)/X)
      EMSQ=DOT(P,5,5)
      Y12=2*DOT(P,ABS(I),ABS(J))/EMSQ*SIGN(1,I*J)
      Y13=2*DOT(P,ABS(I),ABS(K))/EMSQ*SIGN(1,I*K)
      Y23=2*DOT(P,ABS(J),ABS(K))/EMSQ*SIGN(1,J*K)
      L12=LOG(ABS(Y12))
      L13=LOG(ABS(Y13))
      L23=LOG(ABS(Y23))
C---FIRST DECOMPOSE Q AS A*P1+B*P2+C*P3+D*E_(MU,P1,P2,P3)/EMSQ
      P1Q=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)-P(2,ABS(I))*Q(2)
     $     -P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2Q=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)-P(2,ABS(J))*Q(2)
     $     -P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3Q=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)-P(2,ABS(K))*Q(2)
     $     -P(1,ABS(K))*Q(1))*SIGN(1,K)
      QQ=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      A=(P3Q*Y12+P2Q*Y13-P1Q*Y23)/(Y12*Y13*EMSQ)
      B=(P1Q*Y23+P3Q*Y12-P2Q*Y13)/(Y12*Y23*EMSQ)
      C=(P2Q*Y13+P1Q*Y23-P3Q*Y12)/(Y23*Y13*EMSQ)
      D=SQRT(MAX(Z,4*(QQ/EMSQ-A*B*Y12-A*C*Y13-B*C*Y23)/(-Y12*Y13*Y23)))
C---THEN CALCULATE THE CONTRACTIONS WITH -P1(MU)P1(NU) ETC
      R1312=R(Y13,Y12)
      R2312=R(Y23,Y12)
      R2313=R(Y23,Y13)
      CF1=R1312*Y12-R2312*Y12**2/Y13
     $     -L12*Y12*Y13/(Y13+Y23)**2-L23*Y12+Y12/2*(Y23-Y13)/(Y23+Y13)
      CA1=-R1312*Y12/2+R2312*Y12**2/Y13/2+R2313*Y12/2
     $     +L12*Y12*Y13/(Y13+Y23)**2/2-Y12*Y23/(Y23+Y13)/2
      CF2=R2312*Y12-R1312*Y12**2/Y23
     $     -L12*Y12*Y23/(Y13+Y23)**2-L13*Y12+Y12/2*(Y13-Y23)/(Y13+Y23)
      CA2=-R2312*Y12/2+R1312*Y12**2/Y23/2+R2313*Y12/2
     $     +L12*Y12*Y23/(Y13+Y23)**2/2-Y12*Y13/(Y13+Y23)/2
      CF3=R1312*Y12+R2312*Y12-2*L12*Y12
     $     -L13*Y12*Y23/2/(Y12+Y23)**2*(1+Y12+Y23)
     $     -L23*Y12*Y13/2/(Y12+Y13)**2*(1+Y12+Y13)
     $     -Y12/2*(Y13/(Y12+Y13)+Y23/(Y12+Y23))
      CA3=-R1312*Y12/2-R2312*Y12/2+R2313*Y12
     $     +L12*Y12+L13*Y12*Y13/(Y12+Y23)/2+L23*Y12*Y23/(Y12+Y13)/2
      CF4=R1312*Y12/4/Y23*(Y12-Y12**2-Y12*Y23+2*Y12**2*Y23
     $                         +2*Y12*Y23**2+Y23**3)
     $     +R2312*Y12/4/Y13*(Y12-Y12**2-Y12*Y13+2*Y12**2*Y13
     $                           +2*Y12*Y13**2+Y13**3)
     $     +L12*Y12**2/4/(Y13+Y23)*(Y13+Y23-2*Y13*Y23)
     $     -L13*Y12/8/(Y12+Y23)*Y13*(Y12+Y13)*(Y23-2*Y12)
     $     -L23*Y12/8/(Y12+Y13)*Y23*(Y12+Y23)*(Y13-2*Y12)
     $     +Y12/8*(Y13+Y23-2*Y13*Y23)
      CA4=-R1312*Y12/8/Y23*(Y12-Y12**2-Y12*Y23+2*Y12**2*Y23
     $                          +2*Y12*Y23**2+Y23**3)
     $     -R2312*Y12/8/Y13*(Y12-Y12**2-Y12*Y13+2*Y12**2*Y13
     $                           +2*Y12*Y13**2+Y13**3)
     $     +R2313*Y12/8*((1-Y13)**2+(1-Y23)**2)
     $     -L12*Y12**2/8/(Y13+Y23)*(Y13+Y23-2*Y13*Y23)
     $     -L13*Y12/8*Y13*(Y12+Y13)-L23*Y12/8*Y23*(Y12+Y23)
     $     -Y12/8*(Y13+Y23-2*Y13*Y23)
C---AND COMBINE THEM WITH THE APPROPRIATE WEIGHTS
      CFSUB=(A-B)*(A-C)*CF1+(B-C)*(B-A)*CF2+(C-A)*(C-B)*CF3+D*D*CF4
      CASUB=(A-B)*(A-C)*CA1+(B-C)*(B-A)*CA2+(C-A)*(C-B)*CA3+D*D*CA4
      TFSUB=0
      LEIV=CF*(CF*CFSUB+CA*CASUB)
      CFSUB=CFSUB/LEIV
      CASUB=CASUB/LEIV
      LEIV=16*PISQ*LEIV
      END
C-----------------------------------------------------------------------
      SUBROUTINE KPFUNS(X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
      IMPLICIT NONE
C---EVALUATE THE SUM OF THE K AND P FUNCTIONS.
C   IF X<0, RETURN THE DELTA-FUNCTION AND `PLUS' SUBTRACTIONS FOR -X
C   KQF=-SUM_I T_I.T_Q/T_I.T_I GAMMA_I/C_Q
C   PQF=-SUM_I T_I.T_Q/T_Q.T_Q*LOG(Q^2/2P_I.P_Q)
C   WHERE Q IS AN EXTERNALLY AGREED RENORMALIZATION POINT
C   AND LIKEWISE KGF AND PGF
      DOUBLE PRECISION X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG,Z,L,S,
     $     DIS,DILOG,D,LM
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      Z=ABS(X)
      S=1
      IF (X.LE.0) S=-1
      L=LOG((1-Z)/Z)
C---THE PLUS DISTRIBUTIONS
      QQ=QQ+S*XJAC*CF*2/(1-Z)*(L-KQF/2)
     $     -S*XJAC*CF*(1+Z**2)/(1-Z)*(LOG(SCALE)+PQF)
      GG=GG+S*XJAC*CA*2/(1-Z)*(L-KGF/2)
     $     -S*XJAC*CA*2/(1-Z)*(LOG(SCALE)+PGF)
      IF (SCHEME.NE.0) THEN
        DIS=S*XJAC*CF*((1+Z**2)/(1-Z)*(L-0.75)+0.25*(9+5*Z))
        QQ=QQ-DIS
        QG=QG+DIS
      ENDIF
      IF (X.LE.0) THEN
C---THE DELTA FUNCTIONS
        D=DILOG(1-XMIN)
        LM=LOG(1-XMIN)
        QQ=QQ-CF*(5-PISQ+KQF+PISQ/3-LM**2-2*D+KQF*LM
     $       +(2*LM+XMIN+XMIN**2/2)*(LOG(SCALE)+PQF))
        GG=GG-CA*(50D0/9-PISQ+KGF+PISQ/3-LM**2-2*D+KGF*LM
     $       +2*LM*(LOG(SCALE)+PGF))+TR*NF*16D0/9
     $       -(11D0/6*CA-2D0/3*NF*TR)*(LOG(SCALE)+PGF)
      ELSE
C---THE SMOOTH FUNCTIONS
        QQ=QQ+XJAC*CF*(-(1+Z)*L+(1-Z))
        GQ=GQ+XJAC*TR*((Z**2+(1-Z)**2)*(L-LOG(SCALE)-PQF)+2*Z*(1-Z))
        QG=QG+XJAC*CF*((1+(1-Z)**2)/Z*(L-LOG(SCALE)-PGF)+Z)
        GG=GG+XJAC*CA*((1-Z)/Z-1+Z*(1-Z))*2*(L-LOG(SCALE)-PGF)
        IF (SCHEME.NE.0) THEN
          DIS=XJAC*TR*((Z**2+(1-Z)**2)*L+8*Z*(1-Z)-1)
          GQ=GQ-DIS
          GG=GG+2*NF*DIS
        ENDIF
      ENDIF
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE KPFUNS_SCL_VAR(X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
      IMPLICIT NONE
C---EVALUATE THE SUM OF THE K AND P FUNCTIONS.
C   IF X<0, RETURN THE DELTA-FUNCTION AND `PLUS' SUBTRACTIONS FOR -X
C   KQF=-SUM_I T_I.T_Q/T_I.T_I GAMMA_I/C_Q
C   PQF=-SUM_I T_I.T_Q/T_Q.T_Q*LOG(Q^2/2P_I.P_Q)
C   WHERE Q IS AN EXTERNALLY AGREED RENORMALIZATION POINT
C   AND LIKEWISE KGF AND PGF
      DOUBLE PRECISION X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ(3),GQ(3),QG(3)
     $     ,GG(3),Z,L,S,DIS,DILOG,D,LM
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      LOGICAL SCALE_VAR
      DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
      DOUBLE PRECISION MUF(3)
      DATA MUF /1d0, 4d0, 0.25d0/
      Z=ABS(X)
      S=1
      IF (X.LE.0) S=-1
      L=LOG((1-Z)/Z)
C---THE PLUS DISTRIBUTIONS
      QQ=QQ+S*XJAC*CF*2/(1-Z)*(L-KQF/2)
     $     -S*XJAC*CF*(1+Z**2)/(1-Z)*(LOG(MUF)+PQF)
      GG=GG+S*XJAC*CA*2/(1-Z)*(L-KGF/2)
     $     -S*XJAC*CA*2/(1-Z)*(LOG(MUF)+PGF)
      IF (SCHEME.NE.0) THEN
        DIS=S*XJAC*CF*((1+Z**2)/(1-Z)*(L-0.75)+0.25*(9+5*Z))
        QQ=QQ-DIS
        QG=QG+DIS
      ENDIF
      IF (X.LE.0) THEN
C---THE DELTA FUNCTIONS
        D=DILOG(1-XMIN)
        LM=LOG(1-XMIN)
        QQ=QQ-CF*(5-PISQ+KQF+PISQ/3-LM**2-2*D+KQF*LM
     $       +(2*LM+XMIN+XMIN**2/2)*PQF)
        QQ=QQ-CF*(2*LM+XMIN+XMIN**2/2)*LOG(MUF)
        GG=GG-CA*(50D0/9-PISQ+KGF+PISQ/3-LM**2-2*D+KGF*LM
     $       +2*LM*(LOG(MUF)+PGF))+TR*NF*16D0/9
     $       -(11D0/6*CA-2D0/3*NF*TR)*PGF
        GG=GG-(11D0/6*CA-2D0/3*NF*TR)*LOG(MUF)
      ELSE
C---THE SMOOTH FUNCTIONS
        QQ=QQ+XJAC*CF*(-(1+Z)*L+(1-Z))
        GQ=GQ+XJAC*TR*((Z**2+(1-Z)**2)*(L-LOG(MUF)-PQF)+2*Z*(1-Z))
        QG=QG+XJAC*CF*((1+(1-Z)**2)/Z*(L-LOG(MUF)-PGF)+Z)
        GG=GG+XJAC*CA*((1-Z)/Z-1+Z*(1-Z))*2*(L-LOG(MUF)-PGF)
        IF (SCHEME.NE.0) THEN
          DIS=XJAC*TR*((Z**2+(1-Z)**2)*L+8*Z*(1-Z)-1)
          GQ=GQ-DIS
          GG=GG+2*NF*DIS
        ENDIF
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE COLTHR(S,P,W,*)
      IMPLICIT NONE
C---GENERATE A COLLINEAR SPLITTING TO GIVE THREE PARTONS
C   AND EVALUATE THE WEIGHT FOR IT
      INTEGER I
      DOUBLE PRECISION S,P(4,7),W(-6:6),M(-6:6),O,X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,DOT
      PARAMETER (O=0)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      DOUBLE PRECISION QQscl(3) ,GQscl(3),QGscl(3) ,GGscl(3)
      LOGICAL SCALE_VAR
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTWO(P,M)
C---IN FACT THE GENERATION WAS ALREADY DONE EARLIER
      CALL GETCOL(X,XJAC)
C---AND THE KINEMATICS
      DO I=1,4
        P(I,1)=P(I,1)/X
        P(I,3)=P(I,1)*(1-X)
      ENDDO
      return !AK: Just need the kinematics 
C---SO WE JUST HAVE TO CALCULATE THE WEIGHT
      XMIN=2*DOT(P,1,6)/S
      QQ=0
      GG=0
      GQ=0
      QG=0
      KQF=1.5
      if(SCALE_VAR) then
         QQscl = QQ
         GQscl = GQ
         QGscl = QG
         GGscl = GG
         CALL KPFUNS_SCL_VAR(-X,XJAC,XMIN,KQF,O,O,O,QQscl,GQscl,QGscl
     $        ,GGscl)
C---  THE TOTAL
         SCL_WEIGHT(:,0)=GGscl(:)*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               SCL_WEIGHT(:,I)=QGscl(:)*M(0)+QQscl(:) *M(I)
               IF (ABS(I).LE.NF) SCL_WEIGHT(:,0)=GQscl(:)*M(I)
            ENDIF
         ENDDO
         W(:) = SCL_WEIGHT(1,:) 
         do i = 1,3
            SCL_WEIGHT(i,:) = SCL_WEIGHT(i,:) / W(:) 
         enddo
      else
         CALL KPFUNS(-X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---  THE TOTAL
         W(0)=GG*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               W(I)=QG*M(0)+QQ*M(I)
               IF (ABS(I).LE.NF) W(0)=GQ*M(I)
            ENDIF
         ENDDO
      endif
!      CALL KPFUNS(X,XJAC,XMIN,KQF,O,O,O,QQ,GQ,QG,GG)
C---THE TOTAL
!      W(0)=GG*M(0)
!      DO I=-6,6
!        IF (I.NE.0) THEN
!          W(I)=QG*M(0)+QQ*M(I)
!          IF (ABS(I).LE.NF) W(0)=W(0)+GQ*M(I)
!        ENDIF
!      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE COLFOR(S,P,W,*)
      IMPLICIT NONE
C---GENERATE A COLLINEAR SPLITTING TO GIVE FOUR PARTONS
C   AND EVALUATE THE WEIGHT FOR IT
      INTEGER I, INF
      DOUBLE PRECISION S,P(4,7),W(-6:6),M(-6:6),X,XJAC,XMIN,
     $     QQ,GQ,QG,GG,KQF,KGF,PQF,PGF,L12,L13,DOT,EMSQ
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
      DOUBLE PRECISION QQscl(3) ,GQscl(3),QGscl(3) ,GGscl(3)
      LOGICAL SCALE_VAR
      COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR
C---CALCULATE THE LOWEST-ORDER MATRIX-ELEMENT
      CALL MATTHR(P,M)
C---IN FACT THE GENERATION WAS ALREADY DONE EARLIER
      CALL GETCOL(X,XJAC)
C---SO WE JUST HAVE TO CALCULATE THE WEIGHT
      XMIN=2*DOT(P,1,6)/S
      QQ=0
      QG=0
      GQ=0
      GG=0
      KQF=(1.5*(CF-CA/2)+0.5*(11D0/6*CA-2D0/3*NF*TR))/CF
      KGF=1.5
      EMSQ=-DOT(P,5,5)
      L12=LOG(2*DOT(P,1,2)/EMSQ)
      L13=LOG(2*DOT(P,1,3)/EMSQ)
      PQF=-((CF-CA/2)/CF*L12+CA/2/CF*L13)
      PGF=-(L12+L13)/2
      if(SCALE_VAR) then
         QQscl = QQ
         GQscl = GQ
         QGscl = QG
         GGscl = GG
         CALL KPFUNS_SCL_VAR(-X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQscl,GQscl
     $        ,QGscl,GGscl)
C---  THE TOTAL
         SCL_WEIGHT(:,0)=GGscl(:)*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               SCL_WEIGHT(:,I)=QGscl(:)*M(0)+QQscl(:) *M(I)
               IF (ABS(I).LE.NF) SCL_WEIGHT(:,0)=GQscl(:)*M(I)
            ENDIF
         ENDDO
         ! AK: Something is fishy here... Why is W(0) = 0?? 
         W(:) = SCL_WEIGHT(1,:) 
         do i = 1,3
            do inf = -6,6
               if(W(inf) .ne.0d0) then
                  SCL_WEIGHT(i,inf) = SCL_WEIGHT(i,inf) / W(inf)
               endif
            enddo
         enddo
      else
         CALL KPFUNS(-X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
C---  THE TOTAL
         W(0)=GG*M(0)
         DO I=-6,6
            IF (I.NE.0) THEN
               W(I)=QG*M(0)+QQ*M(I)
               IF (ABS(I).LE.NF) W(0)=GQ*M(I)
            ENDIF
         ENDDO
      endif
!      CALL KPFUNS(X,XJAC,XMIN,KQF,KGF,PQF,PGF,QQ,GQ,QG,GG)
C---THE TOTAL
!      W(0)=GG*M(0)
!      DO I=-6,6
!        IF (I.NE.0) THEN
!          W(I)=QG*M(0)+QQ*M(I)
!          IF (ABS(I).LE.NF) W(0)=W(0)+GQ*M(I)
!        ENDIF
!      ENDDO
C---AND THE KINEMATICS
      DO I=1,4
        P(I,1)=P(I,1)/X
        P(I,4)=P(I,1)*(1-X)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE SUBTHR(PERM,SS,P,Q,S,JAC,*)
      IMPLICIT NONE
C---GENERATE A TWO-PARTON STATE FROM A THREE-PARTON STATE,
C   CALCULATE THE JACOBIAN FACTOR FOR THE CORRESPONDING CHANNEL,
C   AND (IF PERM.GT.0) THE APPROXIMATE MATRIX-ELEMENT.
      INTEGER NPERM,PERM,IPERM,IJF,KF,I,J,K,M
      DOUBLE PRECISION SS,P(4,7),Q(4,7),S(-6:6),JAC,
     $     X,Z,DEN,F,QQ,GQ,EMSQ,DOT,XMIN
      PARAMETER (NPERM=1)
      DIMENSION IPERM(3,NPERM),IJF(NPERM),KF(NPERM)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      DATA IPERM/2,3,1/
      DATA IJF/2/
      DATA  KF/1/
      IF (ABS(PERM).GT.NPERM.OR.PERM.EQ.0)
     $     STOP 'PERM TOO BIG IN SUBTHR!'
C---FIND WHICH PARTONS TO TREAT
      I=IPERM(1,ABS(PERM))
      J=IPERM(2,ABS(PERM))
      K=IPERM(3,ABS(PERM))
C---FIND THE KINEMATIC VARIABLES
      EMSQ=2*(-DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
      DEN=2/EMSQ
      X=1/(1+DOT(P,I,J)*DEN)
      Z=DOT(P,I,K)*X*DEN
C--- GPS -- avoid crashes later on
      if (z .lt. cutoff .or. (1d0-z) .lt. cutoff 
     $     .or. (1d0-x) .lt. cutoff .or. x .lt. cutoff) return 1
C---COPY INTO Q, REPLACING K BY KTILDE AND I BY IJTILDE
      DO M=1,4
        Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-(1-X)*P(M,K)
        Q(M,KF(ABS(PERM)))=X*P(M,K)
        Q(M,5)=P(M,5)
        Q(M,6)=P(M,6)
        Q(M,7)=P(M,7)
      ENDDO
C---CALCULATE THE CORRESPONDING JACOBIAN FACTOR
      XMIN=2*DOT(Q,1,6)/SS
C--- GPS to avoid crashes -------
      if (xmin .lt. cutoff) return 1
      JAC=1/(2*NPOW(1)*(Z*(1-Z))**XPOW(1))*(Z**XPOW(1)+(1-Z)**XPOW(1))
     $     *(0.5/(-X*LOG(XMIN))
     $     +0.5*((1-XMIN)/(1-X))**XPOW(1)/(NPOW(1)*(1-XMIN)))
      RETURN ! AK: We only need the jacobian
      IF (PERM.LT.0) RETURN
C---OVERALL FACTOR
      F=16*PISQ/EMSQ
C---CALCULATE WEIGHT FOR FINAL-STATE SPLITTING FUNCTION
      QQ=F*CF*(2/(2-Z-X)-(1+Z))/(1-X)
C---CALCULATE WEIGHTS FOR INITIAL-STATE SPLITTING FUNCTIONS
      QQ=QQ+F*CF*(2/(2-Z-X)-(1+X))/(1-Z)
      GQ=F*TR*(X**2+(1-X)**2)/(Z*(1-Z))
C---MULTIPLY WITH THE LOWEST-ORDER MATRIX ELEMENT
      CALL MATTWO(Q,S)
      S(0)=0
      DO I=1,NF
        S(0)=S(0)+GQ*S(I)
      ENDDO
      DO I=-6,6
        IF (I.NE.0) S(I)=QQ*S(I)
      ENDDO
C---READJUST THE MOMENTA A BIT
      DO M=1,4
        Q(M,3)=P(M,1)-Q(M,1)
        Q(M,1)=P(M,1)
      ENDDO
 999  END
C-----------------------------------------------------------------------
      SUBROUTINE SUBFOR(PERM,SS,P,Q,S,JAC,*)
      IMPLICIT NONE
C---GENERATE A THREE-PARTON STATE FROM A FOUR-PARTON STATE,
C   CALCULATE THE JACOBIAN FACTOR FOR THE CORRESPONDING CHANNEL,
C   AND (IF PERM.GT.0) THE APPROXIMATE MATRIX-ELEMENT.
      INTEGER NPERM,PERM,IPERM,IJF,KF,LF,NPERM3,I,J,K,L,M,n
      DOUBLE PRECISION SS,P(4,7),Q(4,7),S(-6:6),JAC,
     $     X,Z,DEN,QQ,GQ,EMSQ,DOT,XMIN,JTMP,STMP(-6:6),QTMP(4,7),QUSQ,
     $     Y,ZI,ZJ,OY,ZTI,ZTJ,V(4),VV,S1(-6:6),S2(-6:6),SC,S3(-6:6),CUT
     $     ,s4(-6:6),s5(-6:6),s6(-6:6),s7(-6:6),s8(-6:6),gg,qg
     $     ,s9(-6:6),s10(-6:6),s11(-6:6),temp
      PARAMETER (NPERM=6,NPERM3=1)
      DIMENSION IPERM(4,NPERM),IJF(NPERM),KF(NPERM),LF(NPERM)
      INTEGER SCHEME,NF
      DOUBLE PRECISION CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ(-6:6),SCALE
      COMMON  /COLFAC/ CF,CA,TR,PI,PISQ,HF,CUTOFF,EQ,SCALE,SCHEME,NF
      INTEGER NPOW(2)
      DOUBLE PRECISION XPOW(2)
      COMMON  /SAMPLE/ XPOW,NPOW
      DATA IPERM/2,3,4,1, 2,4,3,1, 2,3,1,4, 2,4,1,3, 3,4,2,1, 3,4,1,2/
      DATA IJF/2,2,2,2,3,3/
      DATA  KF/3,3,1,1,2,1/
      DATA  LF/1,1,3,3,1,2/
      IF (ABS(PERM).GT.NPERM.OR.PERM.EQ.0)
     $     STOP 'PERM TOO BIG IN SUBFOR!'

      do i=-6,6
        s(i)=0
        s1(i)=0
        s2(i)=0
        s3(i)=0
        s4(i)=0
        s5(i)=0
        s6(i)=0
        s7(i)=0
        s8(i)=0
        s9(i)=0
        s10(i)=0
        s11(i)=0
      enddo

C---FIND WHICH PARTONS TO TREAT
      I=IPERM(1,ABS(PERM))
      J=IPERM(2,ABS(PERM))
      K=IPERM(3,ABS(PERM))
      L=IPERM(4,ABS(PERM))
C---FIND THE KINEMATIC VARIABLES
      IF (K.EQ.1) THEN
        EMSQ=2*(-DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
        DEN=2/EMSQ
        !write(0,*) DOT(P,I,J)*DEN
C---GPS avoid division by zero: but should one perhaps be in any case
C   worried if DOT(P,I,J)*DEN<0?
        X = (1+DOT(P,I,J)*DEN)
        if (X .LE. CUTOFF) return 1
        X = 1/X
c$$$        X=1/(1+DOT(P,I,J)*DEN)
        Z=DOT(P,I,K)*X*DEN
C---GPS --------------
        if (Z .LE. CUTOFF .OR. 1-Z .LE. CUTOFF .or.(1-x).le.cutoff) then
           RETURN 1
        end if
      ELSE
        EMSQ=2*(DOT(P,I,J)+DOT(P,J,K)+DOT(P,K,I))
        DEN=2/EMSQ
        Y=DOT(P,I,J)*DEN
        ZI=DOT(P,I,K)*DEN
        ZJ=DOT(P,J,K)*DEN
C---GPS-------------------
        if (y .le. cutoff .or. 1-y .le. cutoff .or.
     $       1-zj .le. cutoff .or. 1-zi .le. cutoff) then
           return 1
        end if
        OY=1/(1-Y)
        ZTI=ZI*OY
        ZTJ=ZJ*OY
      ENDIF
C---COPY INTO Q, REPLACING K BY KTILDE AND I BY IJTILDE
      DO M=1,4
        IF (K.EQ.1) THEN
          Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-(1-X)*P(M,K)
          Q(M,KF(ABS(PERM)))=X*P(M,K)
        ELSE
          Q(M,IJF(ABS(PERM)))=P(M,I)+P(M,J)-Y*OY*P(M,K)
          Q(M,KF(ABS(PERM)))=OY*P(M,K)
        ENDIF
        Q(M,LF(ABS(PERM)))=P(M,L)
        Q(M,5)=P(M,5)
        Q(M,6)=P(M,6)
        Q(M,7)=P(M,7)
      ENDDO
C---REENFORCE THE MINIMUM CUTOFF ON ALL PAIR MASSES
      CUT=CUTOFF*(DOT(Q,1,2)+DOT(Q,1,3))
      IF (DOT(Q,1,2).LT.CUT.OR.DOT(Q,1,3).LT.CUT.OR.DOT(Q,2,3).LT.CUT)
     $     RETURN 1
C---CALCULATE THE CORRESPONDING JACOBIAN FACTOR
      JAC=0
      DO M=1,NPERM3
        CALL SUBTHR(-M,SS,Q,QTMP,STMP,JTMP,*999)
        JAC=JAC+JTMP
      ENDDO
      IF (K.EQ.1) THEN
        XMIN=2*DOT(Q,1,6)/SS
        JAC=JAC/(2*NPOW(2)*(Z*(1-Z))**XPOW(2))
     $       *(Z**XPOW(2)+(1-Z)**XPOW(2))
     $       *(0.5/(-X*LOG(XMIN))
     $       +0.5*((1-XMIN)/(1-X))**XPOW(2)/(NPOW(2)*(1-XMIN)))
      ELSE
        IF (ABS(PERM).LE.4) THEN
          JAC=JAC/(NPOW(2)**2*(Y*(1-ZI))**XPOW(2))
        ELSE
          JAC=JAC/(2*NPOW(2)**2*(Y*(1-ZI)*(1-ZJ))**XPOW(2)
     $         /((1-ZI)**XPOW(2)+(1-ZJ)**XPOW(2)))
        ENDIF
        JAC=JAC*2
      ENDIF
      QUSQ=-DOT(P,5,5)
      JAC=JAC*QUSQ/EMSQ
C---INCLUDE A PRIORI CHANNEL WEIGHTS
      JAC=JAC/8
      IF (ABS(PERM).GT.4) JAC=JAC*2
      IF (PERM.LT.0) RETURN
C---CALCULATE WEIGHT FOR QUARK-GLUON SPLITTING FUNCTION
      IF (PERM.LE.4) THEN
        CALL MATTHR(Q,S1)
        IF (K.EQ.1) THEN
          QQ=16*PISQ/EMSQ*(X**2+Z**2)/((1-X)*(1-Z))
        ELSE
          QQ=16*PISQ/EMSQ*(2/(1-ZTI*(1-Y))-(1+ZTI))/Y
        ENDIF
        IF (KF(PERM).EQ.3) THEN
          QQ=QQ*HF*CA
        ELSE
          QQ=QQ*(CF-HF*CA)
        ENDIF
C---SYMMETRY FACTORS
        QQ=QQ/2
        S1(0)=0
        DO M=-6,6
          IF (M.NE.0) S1(M)=QQ*S1(M)
        ENDDO
        DO M=-6,6
          S2(M)=0
          S3(M)=0
        ENDDO
      ELSE
C---CALCULATE WEIGHT FOR GLUON-GLUON SPLITTING FUNCTION
        DO M=-6,6
          S1(M)=0
        ENDDO
        CALL MATTHR(Q,S2)
        IF (K.EQ.1) THEN
          DO M=1,4
            V(M)=Z*P(M,I)-(1-Z)*P(M,J)
          ENDDO
          VV = -2*Z*(1-Z)*DOT(P,I,J)
          CALL CONTHR(Q,V,VV,2,-1,3,SC)
          do m=-6,6
            s3(m)=16*pisq/emsq*eq(m)**2*SC*4*z*(1-z)*hf*ca/2/(1-x)
          enddo
          QQ=16*PISQ/EMSQ*
     $         ((2/(2-Z-X)+2/(Z+1-X)-4)/(1-X)
     $         +(2/(2-Z-X)-(1+X))/(1-Z)+(2/(Z+1-X)-(1+X))/Z)
          QQ=QQ*HF*CA
          QQ=QQ/2
        ELSE
          DO M=1,4
            V(M)=ZTI*P(M,I)-ZTJ*P(M,J)
          ENDDO
          VV = -2*ZTI*ZTJ*DOT(P,I,J)
          CALL CONTHR(Q,V,VV,2,-1,3,SC)
          do m=-6,6
            s3(m)=16*pisq/emsq*eq(m)**2*SC*4*zti*ztj*hf*ca/2/y
          enddo
          QQ=16*PISQ/EMSQ*(2/(1-ZTI*(1-Y))+2/(1-ZTJ*(1-Y))-4)/Y
          QQ=QQ*HF*CA
          QQ=QQ/2
        ENDIF
        S2(0)=0
        DO M=-6,6
          IF (M.NE.0) S2(M)=QQ*S2(M)
        ENDDO
      ENDIF
c---calculate weight for quark-gluon splitting function (all fsr)
      if (perm.ne.1.and.perm.ne.3) then
        call matthr(q,s4)
        if (k.eq.1) then
          gg=16*pisq/emsq*(2/(2-z-x)-(1+z))/(1-x)
          gg=gg*hf*ca
        else
          gg=16*pisq/emsq*(2/(1-zti*(1-y))-(1+zti))/y
          gg=gg*(cf-hf*ca)
        endif
        s4(0)=s4(0)*gg
        do m=-6,6
          if (m.ne.0) s4(m)=0
        enddo
      else
        do m=-6,6
          s4(m)=0
        enddo
      endif
c---calculate weight for gluon-quark splitting function (isr)
      if (perm.eq.3.or.perm.eq.4.or.perm.eq.6) then
         if (perm.eq.3) then
            call matthr(q,s5)
            gq=16*pisq/emsq*(x**2+(1-x)**2)/z/(1-z)
            gq=gq*(cf-hf*ca)/cf*tr
         else
!     BUG-FIX: A bug in DISENT was discovered in arXiv:hep-ph/9912488,
!     but not solved. In was since solved in 2005.10705 and some details
!     were given in 2010.07354.
            if(perm.eq.4) then
               do m=1,4
                  temp=q(m,2)
                  q(m,2)=q(m,3)
                  q(m,3)=temp
               enddo
               call matthr(q,s5)
               do m=1,4
                  temp=q(m,2)
                  q(m,2)=q(m,3)
                  q(m,3)=temp
               enddo
            else
               call matthr(q,s5)
            endif
            gq=16*pisq/emsq*(x**2+(1-x)**2)/z
            gq=gq*hf*ca/cf*tr
         endif
         s5(0)=0
         do m=1,nf
            s5(0)=s5(0)+gq*s5(m)
         enddo
         do m=-6,6
            if (m.ne.0) s5(m)=0
         enddo
      else
         do m=-6,6
            s5(m)=0
         enddo
      endif
c---  calculate weight for gluon-gluon splitting function (isr)
      if (perm.eq.4.or.perm.eq.6) then
         call matthr(q,s6)
        do m=1,4
          v(m)=p(m,i)/z-p(m,j)/(1-z)
        enddo
        VV = -2*DOT(P,I,J)/(Z*(1-Z))
        call conthr(q,v,VV,2,3,-1,sc)
        s7(0)=0
        do m=1,nf
          s7(0)=s7(0)-16*pisq/emsq*eq(m)**2*sc*4*(1-x)/x*hf*ca/(1-z)
     $         *tr/cf
        enddo
        do m=-6,6
          if (m.ne.0) s7(m)=0
        enddo
        gg=16*pisq/emsq*
     $       (2/(2-x-z)-2+2*x*(1-x))/(1-z)
        gg=gg*hf*ca
        s6(0)=s6(0)*gg
        do m=-6,6
          if (m.ne.0) s6(m)=0
        enddo
      else
        do m=-6,6
          s6(m)=0
          s7(m)=0
        enddo
      endif
c---calculate weight for quark-antiquark splitting function (fsr)
      if (perm.ge.5) then
        call matthr(q,s8)
        if (k.eq.1) then
          do m=1,4
            v(m)=z*p(m,i)-(1-z)*p(m,j)
          enddo
          VV = -2*Z*(1-Z)*DOT(P,I,J)
          call conthr(q,v,VV,2,-1,3,sc)
          do m=-6,6
            s9(m)=-16*pisq/emsq*eq(m)**2*sc*4*z*(1-z)
     $           *hf*tr*nf/(1-x)
          enddo
          qq=16*pisq/emsq/(1-x)
          qq=qq*hf*tr*nf
        else
          do m=1,4
            v(m)=zti*p(m,i)-ztj*p(m,j)
          enddo
          VV = -2*ZTI*ZTJ*DOT(P,I,J)
          call conthr(q,v,VV,2,-1,3,sc)
          do m=-6,6
            s9(m)=-16*pisq/emsq*eq(m)**2*sc*4*zti*ztj*hf*tr*nf/y
          enddo
          qq=16*pisq/emsq/y
          qq=qq*hf*tr*nf
        endif
        s8(0)=0
        do m=-6,6
          if (m.ne.0) s8(m)=qq*s8(m)
        enddo
      endif
c---calculate weight for quark-antiquark splitting function (isr)
      if (perm.eq.3.or.perm.eq.4) then
        call matthr(q,s10)
        do m=1,4
          v(m)=p(m,i)/z-p(m,j)/(1-z)
        enddo
        VV = -2*DOT(P,I,J)/(Z*(1-Z))
        call conthr(q,v,VV,2,3,-1,sc)
        do m=-6,6
          s11(m)=0
          do n=1,nf
            s11(m)=s11(m)-
     $           16*pisq/emsq*eq(n)**2*sc*4*(1-x)/x*hf*cf/z
     $           *tr/cf
          enddo
        enddo
        s11(0)=0
        qg=16*pisq/emsq*x/z
        qg=qg*hf*cf
        do m=-6,6
          if (m.ne.0) s10(m)=s10(0)*qg
        enddo
        s10(0)=0
      endif
C---ADD THEM TOGETHER
      DO M=-6,6
        S(M)=S1(M)+S2(M)+S3(M)+s4(m)+s5(m)+s6(m)+s7(m)
     $       +s8(m)+s9(m)+s10(m)+s11(m)
      ENDDO
C---READJUST THE MOMENTA A BIT
      DO M=1,4
        Q(M,4)=P(M,1)-Q(M,1)
        Q(M,1)=P(M,1)
      ENDDO
 999  END
C-----------------------------------------------------------------------
      FUNCTION ERTA(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT A FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTA,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234,S
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      S=S12+S13+S14+S23+S24+S34
      ERTA=(S12*S34**2-S13*S24*S34+S14*S23*S34+3*S12*S23*S34+
     $     3*S12*S14*S34+4*S12**2*S34-S13*S23*S24+2*S12*S23*S24-
     $     S13*S14*S24-2*S12*S13*S24+2*S12**2*S24+S14*S23**2+
     $     2*S12*S23**2+S14**2*S23+4*S12*S14*S23+4*S12**2*S23+
     $     2*S12*S14**2+2*S12*S13*S14+4*S12**2*S14+2*S12**2*S13+
     $     2*S12**3)/(2*S13*S134*S234*S24)+
     $     (S24*S34+S12*S34+S13*S24-S14*S23+S12*S13)/(S13*S134**2)+
     $     2*S23*(S-S13)/(S13*S134*S24)+
     $     S34/(2*S13*S24)
      END
C-----------------------------------------------------------------------
      FUNCTION ERTB(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT B FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTB,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S124,S134,S234,S
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S124=S12+S14+S24
      S134=S13+S14+S34
      S234=S23+S24+S34
      S=S12+S13+S14+S23+S24+S34
      ERTB=(S12*S24*S34+S12*S14*S34-S13*S24**2+S13*S14*S24+
     $     2*S12*S14*S24)/(S13*S134*S23*S14)+
     $     S12*(S+S34)*S124/(S134*S234*S14*S24)-
     $     (2*S13*S24+S14**2+S13*S23+2*S12*S13)/(S13*S134*S14)+
     $     S12*S123*S124/(2*S13*S14*S23*S24)
      END
C-----------------------------------------------------------------------
      FUNCTION ERTC(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT C FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTC,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S134,S234
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S134=S13+S14+S34
      S234=S23+S24+S34
      ERTC=-(5*S12*S34**2+2*S12*S24*S34+2*S12*S23*S34+2*S12*S14*S34+
     $     2*S12*S13*S34+4*S12**2*S34-S13*S24**2+S14*S23*S24+
     $     S13*S23*S24+S13*S14*S24-S12*S14*S24-S13**2*S24-
     $     3*S12*S13*S24-S14*S23**2-S14**2*S23+S13*S14*S23-
     $     3*S12*S14*S23-S12*S13*S23)/(4*S134*S234*S34**2)+
     $     (3*S12*S34**2-3*S13*S24*S34+3*S12*S24*S34+3*S14*S23*S34-
     $     S13*S24**2-S12*S23*S34+6*S12*S14*S34+2*S12*S13*S34-
     $     2*S12**2*S34+S14*S23*S24-3*S13*S23*S24-2*S13*S14*S24+
     $     4*S12*S14*S24+2*S12*S13*S24+3*S14*S23**2+2*S14**2*S23+
     $     2*S14**2*S12+2*S12**2*S14+6*S12*S14*S23-2*S12*S13**2-
     $     2*S12**2*S13)/(4*S13*S134*S234*S34)+
     $     (2*S12*S34**2-2*S13*S24*S34+S12*S24*S34+4*S13*S23*S34+
     $     4*S12*S14*S34+2*S12*S13*S34+2*S12**2*S34-S13*S24**2+
     $     3*S14*S23*S24+4*S13*S23*S24-2*S13*S14*S24+4*S12*S14*S24+
     $     2*S12*S13*S24+2*S14*S23**2+4*S13*S23**2+2*S13*S14*S23+
     $     2*S12*S14*S23+4*S12*S13*S23+2*S12*S14**2+4*S12**2*S13+
     $     4*S12*S13*S14+2*S12**2*S14)/(4*S13*S134*S24*S34)
      ERTC=ERTC-
     $     (S12*S34**2-2*S14*S24*S34-2*S13*S24*S34-S14*S23*S34+
     $     S13*S23*S34+S12*S14*S34+2*S12*S13*S34-2*S14**2*S24-
     $     4*S13*S14*S24-4*S13**2*S24-S14**2*S23-S13**2*S23+
     $     S12*S13*S14-S12*S13**2)/(2*S13*S34*S134**2)+
     $     (S12*S34**2-4*S14*S24*S34-2*S13*S24*S34-2*S14*S23*S34-
     $     4*S13*S23*S34-4*S12*S14*S34-4*S12*S13*S34-2*S13*S14*S24+
     $     2*S13**2*S24+2*S14**2*S23-2*S13*S14*S23-S12*S14**2-
     $     6*S12*S13*S14-S12*S13**2)/(4*S34**2*S134**2)
      END
C-----------------------------------------------------------------------
      FUNCTION ERTD(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT D FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTD,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S134
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S134=S13+S14+S34
      ERTD=(S13*S23*S34+S12*S23*S34-S12**2*S34+S13*S23*S24+
     $     2*S12*S23*S24-S14*S23**2+S12*S13*S24+S12*S14*S23+
     $     S12*S13*S14)/(S13**2*S123**2)-
     $     (S12*S34**2-S13*S24*S34+S12*S24*S34-S14*S23*S34-
     $     S12*S23*S34-S13*S24**2+S14*S23*S24-S13*S23*S24-
     $     S13**2*S24+S14*S23**2)/(S13**2*S123*S134)
      END
C-----------------------------------------------------------------------
      FUNCTION ERTE(P,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE ERT E FUNCTION WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION ERTE,P(4,7),DOT,
     $     S12,S13,S14,S23,S24,S34,S123,S124,S134,S234
      S12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      S13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      S14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      S23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      S24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      S34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      S123=S12+S13+S23
      S124=S12+S14+S24
      S134=S13+S14+S34
      S234=S23+S24+S34
      ERTE=(S12*S23*S34-S12*S24*S34+S12*S14*S34+S12*S13*S34+S13*S24**2-
     $     S14*S23*S24+S13*S23*S24+S13*S14*S24+S13**2*S24-S14*S23**2-
     $     S14**2*S23-S13*S14*S23)/(S13*S23*S123*S134)-
     $     S12*(S12*S34-S23*S24-S13*S24-S14*S23-
     $     S14*S13)/(S13*S23*S123**2)-
     $     (S14+S13)*(S24+S23)*S34/(S13*S23*S134*S234)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIA(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIA,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIA = 0
      LEIA = LEIA + P1K*P2K * (  - 2*Q2*P12*P25 - Q2*P14*
     +     P25 - Q2*P23*P25 + P12*P13*P25 + P12*P24*P25 + 2*P12*P25*
     +     P34 - P13*P23*P25 - 2*P13*P24*P25 - P13*P25*P34 + 2*P14*P15*
     +     P24 - 2*P14*P23*P25 - P14*P24*P25 - 2*P15*P23*P25 - 2*P15*P24
     +     *P25 - 2*P15*P25*P34 - P24*P25*P34 )
      LEIA = LEIA + P1K*P3K * (  - Q2*P12 + Q2*P24 + P12*
     +     P13 + P14*P24 + 2*P15*P24 - P24*P34 )*P25
      LEIA = LEIA + P1K*P4K * (  - Q2*P12 - Q2*P23 + P12*
     +     P13 + P12*P34 - P14*P23 )*P25
      LEIA = LEIA + P1K**2 * ( Q2*P23 + Q2*P24 - P13*P23
     +     - P13*P24 - P24*P34 )*P25
      LEIA = LEIA + P2K*P3K * (  - Q2*P12*P25 - Q2*P14*
     +     P25 - 2*P12*P15*P25 + P12*P24*P25 + P12*P25*P34 + 2*P14*P15*
     +     P24 - P14*P23*P25 - 4*P15*P23*P25 - 2*P15*P24*P25 - 2*P15*P25
     +     *P34 )
      LEIA = LEIA + P2K*P4K * (  - Q2*P12*P25 + Q2*P13*
     +     P25 + P12*P24*P25 + 2*P13*P15*P25 + P13*P23*P25 - P13*P25*P34
     +     + 2*P14*P15*P24 - 2*P15*P23*P25 - 2*P15*P24*P25 )
      LEIA = LEIA + P2K**2 * ( Q2*P13 + Q2*P14 + 2*P13*P15
     +     - P13*P24 - P13*P34 - P14*P24 + 2*P15*P34 )*P25
      LEIA = LEIA + P3K*P4K * (  - 2*P12*P15 + P12*P34 - P13*
     +     P24 - P14*P23 - 2*P15*P23 - P15*P25 )*P25
      LEIA = LEIA + P3K**2 * ( P14 + 2*P15 )*P24*P25
      LEIA = LEIA + P4K**2 * ( P13*P23*P25 )
c$$$      LEIA = LEIA + KK * ( 2*Q2*P12*P14*P25 + 2*Q2*P12*P23
c$$$     +     *P25 + 2*Q2*P12**2*P25 - 2*Q2*P14*P15*P24 + 2*Q2*P14*
c$$$     +     P23*P25 + P12*P13*P25*P34 + 4*P12*P15*P23*P25 + 2*P12*P15*P25
c$$$     +     *P34 + P12*P24*P25*P34 - P13*P14*P23*P25 - 2*P13*P14*P24*P25
c$$$     +     - 2*P13*P15*P24*P25 - 2*P13*P23*P24*P25 - P13*P24**2*P25 - 
c$$$     +     P13**2*P24*P25 + 2*P14*P15*P23*P25 - P14*P23*P24*P25 + 4*P15*
c$$$     +     P23*P24*P25 + 4*P15*P23*P25*P34 + 4*P15*P23**2*P25 + 2*P15*
c$$$     +     P24*P25*P35 + P15*P25**2*P34 )/4
      LEIA = LEIA/(P15*P25**2*P13*P24)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIB(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIB,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIB = 0
      LEIB = LEIB + P1K*P2K * (  - 4*Q2*P14*P24*P45 - 2*P12
     +     *P15*P25*P34 - P12*P15*P25*P35 - P12*P15*P25*P45 + 8*P13*P15*
     +     P23*P24 + 4*P13*P15*P24*P34 + 4*P13*P15*P24**2 + 4*P14*P15*
     +     P23*P24 + 2*P14*P15*P24*P45 + 2*P14*P24*P25*P45 - 4*P15*P23*
     +     P24*P25 + 4*P15*P23*P24*P34 )
      LEIB = LEIB + P1K*P3K * (  - P12*P14*P15*P25 - 2*P12*
     +     P14*P24*P45 + P12*P15*P24*P25 + 4*P12*P15*P24*P34 - P12**2*
     +     P15*P25 - 4*P13*P15*P24**2 + 4*P14*P15*P23*P24 )
      LEIB = LEIB + P1K*P4K * (  - P12*P13*P15*P25 - 2*P12*
     +     P14*P24*P45 + P12*P15*P23*P25 - P12**2*P15*P25 + 4*P13*P15*
     +     P23*P24 + 4*P14*P15*P23*P24 )
      LEIB = LEIB + P1K**2 * ( 2*Q2*P14*P24*P45 - 2*P12*P14*
     +     P24*P45 + P12*P15**2*P25 - 2*P14*P24*P25*P45 - 4*P15*P23*P24*
     +     P34 )
      LEIB = LEIB + P2K*P3K * ( 4*P12*P13*P15*P24 + P12*P14*
     +     P15*P25 - 2*P12*P14*P24*P45 - P12*P15*P24*P25 - P12**2*P15*
     +     P25 + 4*P13*P15*P23*P24 + 4*P13*P15*P24**2 )
      LEIB = LEIB + P2K*P4K * (  - 4*P12*P13*P15*P24 + P12*
     +     P13*P15*P25 - 8*P12*P14*P15*P24 - 2*P12*P14*P24*P45 - P12*P15
     +     *P23*P25 - 4*P12*P15*P24*P34 - P12**2*P15*P25 - 4*P13*P14*P15
     +     *P24 + 4*P13*P15*P23*P24 + 4*P13*P15*P24**2 - 4*P13**2*P15*
     +     P24 )
      LEIB = LEIB + P2K**2 * ( 2*Q2*P14*P24*P45 - 2*P12*P14*
     +     P24*P45 + P12*P15*P25**2 - 4*P13*P14*P15*P24 - 4*P13*P15*P24*
     +     P34 - 4*P13**2*P15*P24 - 2*P14*P15*P24*P45 )
      LEIB = LEIB + P3K*P4K * (  - 2*P12*P25 - 4*P14*P24 )
     +     *P12*P15
      LEIB = LEIB + P3K**2 * (  - 4*P12*P14*P15*P24 )
c$$$      LEIB = LEIB + KK * ( 2*Q2*P12*P14*P24*P45 + Q2*
c$$$     +     P12**2*P15*P25 - 2*Q2*P13*P15*P23*P24 - 2*Q2*P14*P15*P23*
c$$$     +     P24 + 2*Q2*P15*P23*P24*P25 + P12*P13*P14*P15*P25 - 4*P12*
c$$$     +     P13*P15*P23*P24 - 2*P12*P13*P15*P24*P34 + 2*P12*P14*P15*P24*
c$$$     +     P34 + 4*P12*P14*P15*P24**2 + P12*P15*P23*P24*P25 - 2*P12*P15*
c$$$     +     P23*P24*P34 + 2*P12*P15*P24**2*P34 - 2*P13*P14*P15*P23*P24 + 
c$$$     +     2*P13*P14*P15*P24**2 - 2*P13*P15*P23*P24**2 - 2*P13*P15*
c$$$     +     P24**3 + 2*P13**2*P15*P24**2 + 2*P14*P15*P23*P24**2 + 2*P14*
c$$$     +     P15*P23**2*P24 - 2*P14**2*P15*P23*P24 - 2*P15**2*P23*P24*P25
c$$$     +     )/2
      LEIB = LEIB/(2*P13*P14*P15*P23*P24*P25)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIC(P,Q,I,J,K,L)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L
      DOUBLE PRECISION LEIC,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIC = 0
      LEIC = LEIC + P1K*P2K * (  - 4*P12*P13*P15*P25*P34 + 6*
     +     P12*P13*P24*P25*P34 - 2*P12*P14*P15*P25*P34 - 2*P12*P14*P24*
     +     P25*P34 - 2*P12*P15*P25*P34**2 + 2*P12*P24*P25*P34**2 + 4*P13
     +     *P14*P15*P24*P34 - 2*P13*P14*P15*P25*P34 - 3*P13*P14*P23*P24*
     +     P25 - P13*P14*P24**2*P25 + P13*P14**2*P15*P24 - 2*P13*P15*P23
     +     *P25*P34 + P13*P15*P24*P25*P34 - P13*P15*P24*P34**2 - P13*P15
     +     *P25*P34**2 + P13*P23*P24*P25*P34 - 2*P13*P24**2*P25*P34 + 6*
     +     P13**2*P14*P15*P24 - P13**2*P23*P24*P25 + P13**2*P24*P25*P34
     +     - 3*P13**2*P24**2*P25 + P13**3*P15*P24 - P14*P15*P23*P25*P34
     +     - 2*P14*P15*P24*P25*P34 - 2*P14*P15*P25*P34**2 - 3*P14*P23*
     +     P24*P25*P34 + P14*P24*P25*P34**2 - 2*P14*P24**2*P25*P34 - 
     +     P14**2*P15*P25*P34 - P14**2*P24*P25*P34 + P15*P23*P25*P34**2
     +     + P15*P24*P25*P34**2 + P15*P25*P34**3 + P24*P25*P34**3 - 
     +     P24**2*P25*P34**2 )
      LEIC = LEIC + P1K*P3K * ( 3*P12*P13*P34 - P12*P14*P34
     +     - P12*P34**2 + P13*P14*P24 + 2*P13*P15*P34 - P13**2*P24 - 
     +     P14*P15*P34 + 3*P14*P24*P34 + P15*P34**2 - 3*P24*P34**2 )
     +     *P24*P25
      LEIC = LEIC + P1K*P4K * (  - 2*P12*P13*P15*P34 + 3*P12*
     +     P13*P24*P34 - P12*P14*P15*P34 - P12*P14*P24*P34 - 3*P12*P15*
     +     P34**2 + 3*P12*P24*P34**2 - P13*P14*P23*P24 - P13*P15*P23*P34
     +     - 2*P13*P15*P24*P34 - 2*P13*P23*P24*P34 + P13**2*P23*P24 + 2
     +     *P14*P15*P23*P34 - 3*P14*P23*P24*P34 + P23*P24*P34**2 )*P25
      LEIC = LEIC + P1K**2 * ( 2*P13*P15 - 3*P13*P23 - 3*P13*
     +     P24 + P14*P15 + P14*P23 + P14*P24 + 3*P15*P34 + P23*P34 - 3*
     +     P24*P34 )*P24*P25*P34
      LEIC = LEIC + P2K*P3K * (  - 2*P12*P13*P15*P25*P34 + 
     +     P12*P13*P24*P25*P34 - P12*P14*P15*P25*P34 - 3*P12*P14*P24*P25
     +     *P34 - P12*P15*P25*P34**2 + P12*P24*P25*P34**2 + 2*P13*P14*
     +     P15*P24*P34 - P13*P14*P15*P25*P34 - P13*P14*P23*P24*P25 - 2*
     +     P13*P14*P24*P25*P34 + P13*P14*P24**2*P25 - 2*P13*P14**2*P15*
     +     P24 - 4*P13*P15*P23*P25*P34 - 2*P13*P15*P25*P34**2 + 2*P13**2
     +     *P14*P15*P24 - 2*P14*P15*P23*P25*P34 - 3*P14*P15*P24*P25*P34
     +     - 2*P14*P15*P25*P34**2 - 3*P14*P23*P24*P25*P34 - P14*P24**2*
     +     P25*P34 - P14**2*P24*P25*P34 )
      LEIC = LEIC + P2K*P4K * (  - P12*P13*P15*P25*P34 + 2*
     +     P12*P13*P24*P25*P34 - 2*P12*P14*P15*P25*P34 - 2*P12*P14*P24*
     +     P25*P34 - 2*P12*P24*P25*P34**2 + 4*P13*P14*P15*P24*P34 - P13*
     +     P14*P15*P25*P34 + P13*P14*P24*P25*P34 - 2*P13*P15*P23*P25*P34
     +     - 3*P13*P15*P24*P25*P34 - P13*P15*P25*P34**2 + 3*P13*P23*P24
     +     *P25*P34 - 2*P13*P24*P25*P34**2 + P13*P24**2*P25*P34 + 2*
     +     P13**2*P14*P15*P24 - 2*P13**2*P15*P24*P34 + 2*P13**2*P15*P25*
     +     P34 + P13**2*P23*P24*P25 - P13**2*P24**2*P25 - 2*P13**3*P15*
     +     P24 - 2*P14*P15*P23*P25*P34 - 4*P14*P15*P24*P25*P34 )
      LEIC = LEIC + P2K**2 * ( 2*P13*P14*P15 + P13*P14*P24 + 
     +     P13*P15*P34 - P13*P24*P34 + 2*P13**2*P15 - P13**2*P24 + 2*P14
     +     *P15*P34 + 2*P14*P24*P34 + 2*P14**2*P15 + 2*P14**2*P24 )
     +     *P25*P34
      LEIC = LEIC + P3K*P4K * (  - P12*P13*P15*P34 + 2*P12*
     +     P13*P24*P34 - P12*P14*P15*P34 - P12*P15*P34**2 + 2*P12*P24*
     +     P34**2 - 2*P13*P14*P23*P24 - 2*P13*P15*P23*P34 - 2*P13*P24**2
     +     *P34 - 2*P13**2*P24**2 - 2*P14*P23*P24*P34 )*P25
      LEIC = LEIC + P3K**2 * ( 2*P13*P14*P24 + 2*P13*P15*P34
     +     + 2*P14*P24*P34 )*P24*P25
      LEIC = LEIC + P4K**2 * ( 2*P12*P15*P34 + 2*P13*P23*P24
     +     + 2*P23*P24*P34 )*P13*P25
c$$$      LEIC = LEIC + KK * (  - 4*P12*P13*P14*P15*P24*P34 + 4*
c$$$     +     P12*P13*P14*P15*P25*P34 + 3*P12*P13*P14*P23*P24*P25 - 2*P12*
c$$$     +     P13*P14*P24*P25*P34 + P12*P13*P14*P24**2*P25 - P12*P13*P14**2
c$$$     +     *P15*P24 + 4*P12*P13*P15*P23*P25*P34 + P12*P13*P15*P24*P34**2
c$$$     +     + 2*P12*P13*P15*P25*P34**2 - 2*P12*P13*P23*P24*P25*P34 - 2*
c$$$     +     P12*P13*P24*P25*P34**2 - 6*P12*P13**2*P14*P15*P24 + P12*
c$$$     +     P13**2*P23*P24*P25 - 4*P12*P13**2*P24*P25*P34 + 3*P12*P13**2*
c$$$     +     P24**2*P25 - P12*P13**3*P15*P24 + 2*P12*P14*P15*P23*P25*P34
c$$$     +     + 4*P12*P14*P15*P24*P25*P34 + 5*P12*P14*P15*P25*P34**2 + 6*
c$$$     +     P12*P14*P23*P24*P25*P34 + P12*P14*P24*P25*P34**2 + 4*P12*P14*
c$$$     +     P24**2*P25*P34 + 2*P12*P14**2*P15*P25*P34 + 2*P12*P14**2*P24*
c$$$     +     P25*P34 - P12*P15*P24*P25*P34**2 + P12*P15*P25*P34**3 - P12*
c$$$     +     P23*P24*P25*P34**2 + 3*P12*P24**2*P25*P34**2 + 4*P12**2*P13*
c$$$     +     P15*P25*P34 - 6*P12**2*P13*P24*P25*P34 + 2*P12**2*P14*P15*P25
c$$$     +     *P34 + 2*P12**2*P14*P24*P25*P34 + 2*P12**2*P15*P25*P34**2 - 2
c$$$     +     *P12**2*P24*P25*P34**2 )/4
c$$$      LEIC = LEIC + KK * (  - 2*P13*P14*P15*P23*P24*P34 + 2*
c$$$     +     P13*P14*P15*P23*P25*P34 + P13*P14*P15*P24*P25*P34 - 4*P13*P14
c$$$     +     *P15*P24**2*P34 + 3*P13*P14*P23*P24*P25*P34 - P13*P14*P23*
c$$$     +     P24**2*P25 + P13*P14*P23**2*P24*P25 - 3*P13*P14*P24**2*P25*
c$$$     +     P34 + 2*P13*P14**2*P15*P23*P24 + P13*P14**2*P23*P24*P25 + 2*
c$$$     +     P13*P15*P23*P24*P25*P34 + 4*P13*P15*P23*P25*P34**2 + 4*P13*
c$$$     +     P15*P23**2*P25*P34 - P13*P15*P24*P25*P34**2 + 3*P13*P15*
c$$$     +     P24**2*P25*P34 - 3*P13*P23*P24**2*P25*P34 + 2*P13*P24**2*P25*
c$$$     +     P34**2 - P13*P24**3*P25*P34 - 2*P13**2*P14*P15*P23*P24 - 2*
c$$$     +     P13**2*P14*P15*P24**2 - P13**2*P14*P23*P24*P25 - P13**2*P14*
c$$$     +     P24**2*P25 - 2*P13**2*P15*P24*P25*P34 + 2*P13**2*P15*P24**2*
c$$$     +     P34 - P13**2*P23*P24**2*P25 + P13**2*P24**2*P25*P34 + P13**2*
c$$$     +     P24**3*P25 + 2*P13**3*P15*P24**2 + P13**3*P24**2*P25 + 5*P14*
c$$$     +     P15*P23*P24*P25*P34 + P14*P15*P23*P25*P34**2 + 2*P14*P15*
c$$$     +     P23**2*P25*P34 + 4*P14*P15*P24**2*P25*P34 + 2*P14*P23*P24*P25
c$$$     +     *P34**2 + )/4
c$$$      LEIC = LEIC + KK * ( P14*P23*P24**2*P25*P34 + 3*P14*
c$$$     +     P23**2*P24*P25*P34 - P14**2*P15*P23*P25*P34 + 3*P14**2*P23*
c$$$     +     P24*P25*P34 )/4
      LEIC = LEIC/(2*P15*P25**2*P13*P34**2*P24)
      END
C-----------------------------------------------------------------------
      FUNCTION LEID(P,Q,II,JJ,KKK,LL)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L,II,JJ,KKK,LL
      DOUBLE PRECISION LEID,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
C---CONVERT FROM ERT NOTATION (Q'QQ'BARQBAR) TO LEIDEN (QQBARQ'BARQ')
      I=JJ
      J=LL
      K=KKK
      L=II
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEID = 0
      LEID = LEID + P1K*P2K * (  - 2*P12*P25*P34 - 4*P13*P14*
     +     P15 - 2*P13*P15*P34 + 2*P13*P24*P25 - 2*P14*P15*P34 + 2*P14*
     +     P23*P25 )/P25
      LEID = LEID + P1K*P3K * (  - P12*P34 - 2*P14*P24 + P24*P25 )
      LEID = LEID + P1K*P4K * (  - P12*P34 - 2*P13*P23 + P23*P25 )
      LEID = LEID + P1K**2 * ( P23 + P24 )*P34
      LEID = LEID + P2K*P3K * (  - P12*P25*P34 - 4*P13*P14*P15
     +     - 2*P13*P15*P34 + 3*P14*P15*P25 - 2*P14*P15*P34 - 2*P14*P24*
     +     P25 )/P25
      LEID = LEID + P2K*P4K * (  - P12*P25*P34 - 4*P13*P14*P15
     +     + 3*P13*P15*P25 - 2*P13*P15*P34 - 2*P13*P23*P25 - 2*P14*P15*
     +     P34 )/P25
      LEID = LEID + P2K**2 * ( P13 + P14 )*P34
      LEID = LEID + P3K*P4K * (  - 2*P12*P34 + 2*P13*P24 + 2*
     +     P14*P23 )
      LEID = LEID + P3K**2 * (  - 2*P14*P24 )
      LEID = LEID + P4K**2 * (  - 2*P13*P23 )
c$$$      LEID = LEID + KK * ( Q2*P12*P25*P34 - Q2*P13*P24*P25
c$$$     +     - Q2*P14*P23*P25 + 2*P12*P13*P15*P25 - P12*P13*P24*P25 - 2
c$$$     +     *P12*P13**2*P15 + 2*P12*P14*P15*P25 - P12*P14*P23*P25 - 2*P12
c$$$     +     *P14**2*P15 + P12*P25*P34**2 + P12**2*P25*P34 + 2*P13*P14*P23
c$$$     +     *P25 + 2*P13*P14*P24*P25 + 2*P13*P15*P23*P25 + 2*P13*P23*P24*
c$$$     +     P25 - P13*P24*P25*P34 - 2*P13**2*P15*P23 - 2*P13**2*P15*P24
c$$$     +     + 2*P14*P15*P24*P25 + 2*P14*P23*P24*P25 - P14*P23*P25*P34 - 
c$$$     +     2*P14**2*P15*P23 - 2*P14**2*P15*P24 )/(4*P25)
      LEID = LEID/(P15*P25*P34**2)
      END
C-----------------------------------------------------------------------
      FUNCTION LEIE(P,Q,II,JJ,KKK,LL)
      IMPLICIT NONE
C---EVALUATE THE MATRIX ELEMENT CONTRACTED WITH A TENSOR -Q(MU)Q(NU)
C   WITH P(*,I)=P1, P(*,J)=P2, ETC
      INTEGER I,J,K,L,II,JJ,KKK,LL
      DOUBLE PRECISION LEIE,P(4,7),Q(4),DOT,
     $     P12,P13,P14,P23,P24,P34,P15,P25,P35,P45,Q2,
     $     P1K,P2K,P3K,P4K,KK
C---CONVERT FROM ERT NOTATION (Q'QQ'BARQBAR) TO LEIDEN (QQBARQ'BARQ')
      I=JJ
      J=LL
      K=KKK
      L=II
      P12=2*DOT(P,ABS(I),ABS(J))*SIGN(1,I*J)
      P13=2*DOT(P,ABS(I),ABS(K))*SIGN(1,I*K)
      P14=2*DOT(P,ABS(I),ABS(L))*SIGN(1,I*L)
      P23=2*DOT(P,ABS(J),ABS(K))*SIGN(1,J*K)
      P24=2*DOT(P,ABS(J),ABS(L))*SIGN(1,J*L)
      P34=2*DOT(P,ABS(K),ABS(L))*SIGN(1,K*L)
      P15=P23+P24+P34
      P25=P13+P14+P34
      P35=P12+P14+P24
      P45=P12+P13+P23
      Q2=P12+P13+P14+P23+P24+P34
      P1K=(P(4,ABS(I))*Q(4)-P(3,ABS(I))*Q(3)
     $    -P(2,ABS(I))*Q(2)-P(1,ABS(I))*Q(1))*SIGN(1,I)
      P2K=(P(4,ABS(J))*Q(4)-P(3,ABS(J))*Q(3)
     $    -P(2,ABS(J))*Q(2)-P(1,ABS(J))*Q(1))*SIGN(1,J)
      P3K=(P(4,ABS(K))*Q(4)-P(3,ABS(K))*Q(3)
     $    -P(2,ABS(K))*Q(2)-P(1,ABS(K))*Q(1))*SIGN(1,K)
      P4K=(P(4,ABS(L))*Q(4)-P(3,ABS(L))*Q(3)
     $    -P(2,ABS(L))*Q(2)-P(1,ABS(L))*Q(1))*SIGN(1,L)
      KK=Q(4)*Q(4)-Q(3)*Q(3)-Q(2)*Q(2)-Q(1)*Q(1)
      LEIE = 0
      LEIE = LEIE + P1K*P2K * (  - 2*Q2*P13*P15*P25 + Q2*
     +     P23*P25**2 + 2*P13*P14*P15*P25 - 4*P14*P15*P25*P45 - P14*P23*
     +     P25**2 + 2*P14**2*P15*P45 - 2*P15*P23*P25*P34 + 2*P15*P25**2*
     +     P45 - P23*P25**2*P45 )/(P13*P25*P45)
      LEIE = LEIE + P1K*P3K * (  - 2*Q2*P13*P15 + Q2*P23*
     +     P25 - 2*P12*P14*P15 - P14*P23*P25 - 2*P15*P23*P34 + 2*P15*P25
     +     *P45 - P23*P25*P45 )/(P13*P45)
      LEIE = LEIE + P1K*P4K * ( 2*Q2*P25 + 2*P13*P15 - 2*P14*
     +     P25 - 2*P15*P25 - 2*P15*P34 - 2*P25*P45 )*P23/(P13*P45)
      LEIE = LEIE + P1K**2 * (  - 2*P15*P23*P34 )/(P13*P45)
      LEIE = LEIE + P2K*P3K * (  - 2*P12*P15*P25 + 2*P13*P15*
     +     P25 + 2*P14*P15*P45 - 2*P23*P25**2 )*P14/(P13*P25*P45)
      LEIE = LEIE + P2K*P4K * ( Q2*P23*P25**2 + 2*P12*P15*P25
     +     *P34 - 2*P13*P15*P24*P25 + 2*P14*P15*P23*P25 - 2*P14*P15*P25*
     +     P45 - P14*P23*P25**2 + 2*P14**2*P15*P45 - P15*P23*P25**2 )
     +     /(P13*P25*P45)
      LEIE = LEIE + P2K**2 * ( 2*P13*P15 - P23*P25 )*P14/(P13*P45)
      LEIE = LEIE + P3K*P4K * (  - 2*Q2*P13*P15 + Q2*P23*
     +     P25 - 2*P12*P14*P15 + 2*P13*P15*P23 - P14*P23*P25 - 3*P15*P23
     +     *P25 + 2*P15*P25*P45 )/(P13*P45)
      LEIE = LEIE + P3K**2 * (  - 2*P12*P15 - P23*P25 )*P14/(P13*P45)
      LEIE = LEIE + P4K**2 * ( 2*P15*P23 )/P45
c$$$      LEIE = LEIE + KK * ( Q2*P13*P15*P24*P25 + Q2*P14*P23*
c$$$     +     P25**2 - Q2*P14**2*P15*P45 + Q2*P15*P23*P25**2 + Q2*P23
c$$$     +     *P25**2*P45 - Q2**2*P23*P25**2 - P12*P13*P15*P25*P34 + P12
c$$$     +     *P14*P15*P23*P25 - P12*P15*P24*P25*P34 - P12*P15*P25*P34**2 - 
c$$$     +     P12**2*P15*P25*P34 - P13*P14*P15*P23*P25 - P13*P14*P15*P24*
c$$$     +     P25 - P13*P15*P23*P24*P25 - P14*P15*P23*P24*P25 + P14*P15*P23
c$$$     +     *P25*P34 + P14*P15*P25*P35*P45 - P15*P23*P25**2*P45 )
c$$$     +     /(2*P13*P25*P45)
      LEIE = LEIE/(P15*P25*P34)
      END
C-----------------------------------------------------------------------





C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION DILOG(X)
      IMPLICIT NONE
c$$$C---RETURN THE DILOGARITHM (MAXIMUM ERROR = 0.0058, MAX FRAC ERROR = 1%)
c$$$      DOUBLE PRECISION DILOG,X,PISQO6,OX,XX,L2
c$$$      DATA PISQO6,L2/2*0/
c$$$      IF (PISQO6.EQ.0) PISQO6=(ATAN(1D0)*4)**2/6
c$$$      IF (L2.EQ.0) L2=LOG(2D0)
c$$$      IF (X.LT.-1) THEN
c$$$        XX=1/X
c$$$      ELSE
c$$$        XX=X
c$$$      ENDIF
c$$$      IF (XX.LT.-0.5) THEN
c$$$        OX=1+XX
c$$$        DILOG=-PISQO6/2-L2*LOG(-XX)
c$$$     $       -OX**2/4-5*OX**3/24-OX**4/6-131*OX**5/960
c$$$      ELSEIF (XX.LT.0.5) THEN
c$$$        DILOG=XX+XX**2/4+XX**3/9
c$$$      ELSEIF (XX.LT.1) THEN
c$$$        OX=1-XX
c$$$        DILOG=PISQO6-LOG(OX)*LOG(XX)-OX-OX**2/4-OX**3/9
c$$$      ELSEIF (XX.EQ.1) THEN
c$$$        DILOG=PISQO6
c$$$      ELSE
c$$$        WRITE (*,*) 'DILOG CALLED FOR X=',X
c$$$        DILOG=0
c$$$      ENDIF
c$$$      IF (X.LT.-1) DILOG=-DILOG-PISQO6-LOG(-X)**2/2
      double precision dilog,rsp,x,pisqo6
      data pisqo6/0/
      if (pisqo6.eq.0) pisqo6=(atan(1d0)*4)**2/6
      if (x.lt.1) then
        dilog=rsp(x)
      elseif (x.gt.1) then
c---if dilog is complex, return its real part
        dilog=pisqo6-log(x)*log(x-1)-rsp(1-x)
      else
        dilog=pisqo6
      endif
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE RANGEN(N,R)
      IMPLICIT NONE
C---RANDOM NUMBER GENERATOR
C   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C   RETURNS A VECTOR OF N RANDOM VALUES
C   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
C   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
      DOUBLE PRECISION R(*)
      INTEGER N,I,ISEED(2),K,IZ
      DATA ISEED/12345,678900/
      IF (N.LT.0) WRITE (*,'(I10,A,I10,I11)') -N-1,', ISEED=',ISEED
      IF (N.GT.0) THEN
        DO I=1,N
          K=ISEED(1)/53668
          ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
          IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
          K=ISEED(2)/52774
          ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
          IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
          IZ=ISEED(1)-ISEED(2)
          IF (IZ.LT.1) IZ=IZ+2147483562
          R(I)=DBLE(IZ)*4.656613D-10
        ENDDO
      ELSEIF (N.EQ.0) THEN
        ISEED(1)=NINT(R(1))
        ISEED(2)=NINT(R(2))
      ENDIF
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
!     AK: GPS modified routine below
!      SUBROUTINE RANGEN(N,R)
!      IMPLICIT NONE
!C---RANDOM NUMBER GENERATOR
!C   USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
!C   RETURNS A VECTOR OF N RANDOM VALUES
!C   IF (N.EQ.0) THE FIRST TWO VALUES IN R SET THE SEEDS
!C   IF (N.LT.0) PRINT THE CURRENT VALUES OF THE SEEDS
!      DOUBLE PRECISION R(*)
!      INTEGER N,ISEED(3)
!      if (N.GT.0) then
!         call RM48(R,N)
!      else if (N.LT.0) then
!         call RM48UT(iseed(1),iseed(2),iseed(3))
!         WRITE (6,'(I10,A,I10,I11,I11)') -N-1,', ISEED=',
!     $        iseed(1),iseed(2),iseed(3)
!      else ! N=0
!        IF(NINT(R(1)) .eq. 0) then
!           !-- retrieve seed --
!           call RM48UT(iseed(1),iseed(2),iseed(3))
!           R(1) = ISEED(1)
!           R(2) = ISEED(2)
!           R(3) = ISEED(3)
!        else
!!--   set seed -------
!           ISEED(1)=NINT(R(1))
!           ISEED(2)=NINT(R(2))
!           ISEED(3)=NINT(R(3))
!           call RM48IN(iseed(1),iseed(2),iseed(3))
!        end if
!      end if
!      !DATA ISEED/12345,678900/
!      !DATA ISEED/1277158507, 1826842337/
!      !DATA ISEED/1542788427,  474274578/
!c$$$      IF (N.LT.0) WRITE (*,'(I10,A,I10,I11)') -N-1,', ISEED=',ISEED
!c$$$      IF (N.GT.0) THEN
!c$$$        DO I=1,N
!c$$$          K=ISEED(1)/53668
!c$$$          ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
!c$$$          IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
!c$$$          K=ISEED(2)/52774
!c$$$          ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
!c$$$          IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
!c$$$          IZ=ISEED(1)-ISEED(2)
!c$$$          IF (IZ.LT.1) IZ=IZ+2147483562
!c$$$          R(I)=DBLE(IZ)*4.656613D-10
!c$$$        ENDDO
!c$$$      ELSEIF (N.EQ.0) THEN
!c$$$        IF(NINT(R(1)) .eq. 0) then
!c$$$           !-- retrieve seed --
!c$$$           R(1) = ISEED(1)
!c$$$           R(2) = ISEED(2)
!c$$$        else
!c$$$           !-- set seed -------
!c$$$           ISEED(1)=NINT(R(1))
!c$$$           ISEED(2)=NINT(R(2))
!c$$$        end if
!c$$$      ENDIF
!      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     KTCLUS: written by Mike Seymour, July 1992.
C     Last modified February 1994.
C     Please send comments or suggestions to mike@thep.lu.se
C
C     This is a general-purpose kt clustering package.
C     It can handle ee, ep and pp collisions.
C     It is loosely based on the program of Siggi Bethke.
C
C     The time taken (on a 10MIP machine) is (0.2microsec)*N**3
C     where N is the number of particles.
C     Over 90 percent of this time is used in subroutine KTPMIN, which
C     simply finds the minimum member of a one-dimensional array.
C     It is well worth thinking about optimization: on the SPARCstation
C     a factor of two increase was obtained simply by increasing the
C     optimization level from its default value.
C
C     The approach is to separate the different stages of analysis.
C     KTCLUS does all the clustering and records a merging history.
C     It returns a simple list of the y values at which each merging
C     occured. Then the following routines can be called to give extra
C     information on the most recently analysed event.
C     KTCLUR is identical but includes an R parameter, see below.
C     KTYCUT gives the number of jets at each given YCUT value.
C     KTYSUB gives the number of sub-jets at each given YCUT value.
C     KTBEAM gives same info as KTCLUS but only for merges with the beam
C     KTJOIN gives same info as KTCLUS but for merges of sub-jets.
C     KTRECO reconstructs the jet momenta at a given value of YCUT.
C     It also gives information on which jets at scale YCUT belong to
C     which macro-jets at scale YMAC, for studying sub-jet properties.
C     KTINCL reconstructs the jet momenta according to the inclusive jet
C     definition of Ellis and Soper.
C     KTISUB, KTIJOI and KTIREC are like KTYSUB, KTJOIN and KTRECO,
C     except that they only apply to one inclusive jet at a time,
C     with the pt of that jet automatically used for ECUT.
C     KTWICH gives a list of which particles ended up in which jets.
C     KTWCHS gives the same thing, but only for subjets.
C     Note that the numbering of jets used by these two routines is
C     guaranteed to be the same as that used by KTRECO.
C
C     The collision type and analysis type are indicated by the first
C     argument of KTCLUS. IMODE=<TYPE><ANGLE><MONO><RECOM> where
C     TYPE:  1=>ee, 2=>ep with p in -z direction, 3=>pe, 4=>pp
C     ANGLE: 1=>angular kt def., 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
C            where f()=2(cosh(eta)-cos(phi)) is the QCD emission metric
C     MONO:  1=>derive relative pseudoparticle angles from jets
C            2=>monotonic definitions of relative angles
C     RECOM: 1=>E recombination scheme, 2=>pt scheme, 3=>pt**2 scheme
C
C     There are also abbreviated forms for the most common combinations:
C     IMODE=1 => E scheme in e+e-                              (=1111)
C           2 => E scheme in ep                                (=2111)
C           3 => E scheme in pe                                (=3111)
C           4 => E scheme in pp                                (=4111)
C           5 => covariant E scheme in pp                      (=4211)
C           6 => covariant pt-scheme in pp                     (=4212)
C           7 => covariant monotonic pt**2-scheme in pp        (=4223)
C
C     KTRECO no longer needs to reconstruct the momenta according to the
C     same recombination scheme in which they were clustered. Its first
C     argument gives the scheme, taking the same values as RECOM above.
C
C     Note that unlike previous versions, all variables which hold y
C     values have been named in a consistent way:
C     Y()  is the output scale at which jets were merged,
C     YCUT is the input scale at which jets should be counted, and
C          jet-momenta reconstructed etc,
C     YMAC is the input macro-jet scale, used in determining whether
C          or not each jet is a sub-jet.
C     The original scheme defined in our papers is equivalent to always
C     setting YMAC=1.
C     Whenever a YCUT or YMAC variable is used, it is rounded down
C     infinitesimally, so that for example, setting YCUT=Y(2) refers
C     to the scale where the event is 2-jet, even if rounding errors
C     have shifted its value slightly.
C
C     An R parameter can be used in hadron-hadron collisions by
C     calling KTCLUR instead of KTCLUS.  This is as suggested by
C     Ellis and Soper, but implemented slightly differently,
C     as in M.H. Seymour, LU TP 94/2 (submitted to Nucl. Phys. B.).
C     R**2 multiplies the single Kt everywhere it is used.
C     Calling KTCLUR with R=1 is identical to calling KTCLUS.
C     R plays a similar role to the jet radius in a cone-type algorithm,
C     but is scaled up by about 40% (ie R=0.7 in a cone algorithm is
C     similar to this algorithm with R=1).
C     Note that R.EQ.1 must be used for the e+e- and ep versions,
C     and is strongly recommended for the hadron-hadron version.
C     However, R values smaller than 1 have been found to be useful for
C     certain applications, particularly the mass reconstruction of
C     highly-boosted colour-singlets such as high-pt hadronic Ws,
C     as in M.H. Seymour, LU TP 93/8 (to appear in Z. Phys. C.).
C     Situations in which R<1 is useful are likely to also be those in
C     which the inclusive reconstruction method is more useful.
C
C     Also included is a set of routines for doing Lorentz boosts:
C     KTLBST finds the boost matrix to/from the cm frame of a 4-vector
C     KTRROT finds the rotation matrix from one vector to another
C     KTMMUL multiplies together two matrices
C     KTVMUL multiplies a vector by a matrix
C     KTINVT inverts a transformation matrix (nb NOT a general 4 by 4)
C     KTFRAM boosts a list of vectors between two arbitrary frames
C     KTBREI boosts a list of vectors between the lab and Breit frames
C     KTHADR boosts a list of vectors between the lab and hadronic cmf
C       The last two need the momenta in the +z direction of the lepton
C       and hadron beams, and the 4-momentum of the outgoing lepton.
C
C     The main reference is:
C       S. Catani, Yu.L. Dokshitzer, M.H. Seymour and B.R. Webber,
C         Nucl.Phys.B406(1993)187.
C     The ep version was proposed in:
C       S. Catani, Yu.L. Dokshitzer and B.R. Webber,
C         Phys.Lett.285B(1992)291.
C     The inclusive reconstruction method was proposed in:
C       S.D. Ellis and D.E. Soper,
C         Phys.Rev.D48(1993)3160.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE KTCLUS(IMODE,PP,NN,ECUT,Y,*)
      IMPLICIT NONE
C---DO CLUSTER ANALYSIS OF PARTICLES IN PP
C
C   IMODE   = INPUT  : DESCRIBED ABOVE
C   PP(I,J) = INPUT  : 4-MOMENTUM OF Jth PARTICLE: I=1,4 => PX,PY,PZ,E
C   NN      = INPUT  : NUMBER OF PARTICLES
C   ECUT    = INPUT  : DENOMINATOR OF KT MEASURE. IF ZERO, ETOT IS USED
C   Y(J)    = OUTPUT : VALUE OF Y FOR WHICH EVENT CHANGES FROM BEING
C                        J JET TO J-1 JET
C   LAST ARGUMENT IS LABEL TO JUMP TO IF FOR ANY REASON THE EVENT
C   COULD NOT BE PROCESSED (MOST LIKELY DUE TO TOO MANY PARTICLES)
C
C   NOTE THAT THE MOMENTA ARE DECLARED DOUBLE PRECISION,
C   AND ALL OTHER FLOATING POINT VARIABLES ARE DECLARED DOUBLE PRECISION
C
      INTEGER IMODE,NN
      DOUBLE PRECISION PP(4,*)
      DOUBLE PRECISION ECUT,Y(*),ONE
      ONE=1
      CALL KTCLUR(IMODE,PP,NN,ONE,ECUT,Y,*999)
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      SUBROUTINE KTCLUR(IMODE,PP,NN,R,ECUT,Y,*)
      IMPLICIT NONE
C---DO CLUSTER ANALYSIS OF PARTICLES IN PP
C
C   IMODE   = INPUT  : DESCRIBED ABOVE
C   PP(I,J) = INPUT  : 4-MOMENTUM OF Jth PARTICLE: I=1,4 => PX,PY,PZ,E
C   NN      = INPUT  : NUMBER OF PARTICLES
C   R       = INPUT  : ELLIS AND SOPER'S R PARAMETER, SEE ABOVE.
C   ECUT    = INPUT  : DENOMINATOR OF KT MEASURE. IF ZERO, ETOT IS USED
C   Y(J)    = OUTPUT : VALUE OF Y FOR WHICH EVENT CHANGES FROM BEING
C                        J JET TO J-1 JET
C   LAST ARGUMENT IS LABEL TO JUMP TO IF FOR ANY REASON THE EVENT
C   COULD NOT BE PROCESSED (MOST LIKELY DUE TO TOO MANY PARTICLES)
C
C   NOTE THAT THE MOMENTA ARE DECLARED DOUBLE PRECISION,
C   AND ALL OTHER FLOATING POINT VARIABLES ARE DECLARED DOUBLE PRECISION
C
      INTEGER NMAX,IM,IMODE,TYPE,ANGL,MONO,RECO,N,I,J,NN,
     &     IMIN,JMIN,KMIN,NUM,HIST,INJET,IABBR,NABBR
      PARAMETER (NMAX=512,NABBR=7)
      DOUBLE PRECISION PP(4,*)
      DOUBLE PRECISION R,ECUT,Y(*),P,KT,ETOT,RSQ,KTP,KTS,KTPAIR,KTSING,
     &     KTMIN,ETSQ,KTLAST,KTMAX
      LOGICAL FIRST
      CHARACTER TITLE(4,4)*10
C---KT RECORDS THE KT**2 OF EACH MERGING.
C---KTLAST RECORDS FOR EACH MERGING, THE HIGHEST ECUT**2 FOR WHICH THE
C   RESULT IS NOT MERGED WITH THE BEAM (COULD BE LARGER THAN THE
C   KT**2 AT WHICH IT WAS MERGED IF THE KT VALUES ARE NOT MONOTONIC).
C   THIS MAY SOUND POINTLESS, BUT ITS USEFUL FOR DETERMINING WHETHER
C   SUB-JETS SURVIVED TO SCALE Y=YMAC OR NOT.
C---HIST RECORDS MERGING HISTORY:
C   N=>DELETED TRACK N, M*NMAX+N=>MERGED TRACKS M AND N (M<N).
      COMMON /KTCOMM/ETOT,RSQ,P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX),
     &  KT(NMAX),KTLAST(NMAX),HIST(NMAX),NUM
      DIMENSION INJET(NMAX),IABBR(NABBR)
      DATA FIRST,TITLE,IABBR/.TRUE.,
     &     'e+e-      ','ep        ','pe        ','pp        ',
     &     'angle     ','DeltaR    ','f(DeltaR) ','**********',
     &     'no        ','yes       ','**********','**********',
     &     'E         ','Pt        ','Pt**2     ','**********',
     &     1111,2111,3111,4111,4211,4212,4223/
C---CHECK INPUT
      IM=IMODE
      IF (IM.GE.1.AND.IM.LE.NABBR) IM=IABBR(IM)
      TYPE=MOD(IM/1000,10)
      ANGL=MOD(IM/100 ,10)
      MONO=MOD(IM/10  ,10)
      RECO=MOD(IM     ,10)
      IF (NN.GT.NMAX.OR.NN.LT.1.OR.(NN.LT.2.AND.TYPE.EQ.1))
     &     CALL KTWARN('KTCLUS',100,*999)
      IF (TYPE.LT.1.OR.TYPE.GT.4.OR.ANGL.LT.1.OR.ANGL.GT.3.OR.
     &    MONO.LT.1.OR.MONO.GT.2.OR.RECO.LT.1.OR.RECO.GT.3)
     &     CALL KTWARN('KTCLUS',101,*999)
      IF (FIRST) THEN
         WRITE (6,'(/,1X,54(''*'')/A)')
     &   ' KTCLUS: written by Mike Seymour, July 1992.'
         WRITE (6,'(A)')
     &   ' Last modified February 1994.'
         WRITE (6,'(A)')
     &   ' Please send comments or suggestions to mike@thep.lu.se'
         WRITE (6,'(/A,I2,2A)')
     &   '       Collision type =',TYPE,' = ',TITLE(TYPE,1)
         WRITE (6,'(A,I2,2A)')
     &   '     Angular variable =',ANGL,' = ',TITLE(ANGL,2)
         WRITE (6,'(A,I2,2A)')
     &   ' Monotonic definition =',MONO,' = ',TITLE(MONO,3)
         WRITE (6,'(A,I2,2A)')
     &   ' Recombination scheme =',RECO,' = ',TITLE(RECO,4)
         IF (R.NE.1) THEN
         WRITE (6,'(A,F5.2)')
     &   '     Radius parameter =',R
         IF (TYPE.NE.4) WRITE (6,'(A)')
     &   ' R.NE.1 is strongly discouraged for this collision type!'
         ENDIF
         WRITE (6,'(1X,54(''*'')/)')
         FIRST=.FALSE.
      ENDIF
C---COPY PP TO P
      N=NN
      NUM=NN
      CALL KTCOPY(PP,N,P,(RECO.NE.1))
      ETOT=0
      DO 100 I=1,N
         ETOT=ETOT+P(4,I)
 100  CONTINUE
      IF (ETOT.EQ.0) CALL KTWARN('KTCLUS',102,*999)
      IF (ECUT.EQ.0) THEN
         ETSQ=1/ETOT**2
      ELSE
         ETSQ=1/ECUT**2
      ENDIF
      RSQ=R**2
C---CALCULATE ALL PAIR KT's
      DO 210 I=1,N-1
         DO 200 J=I+1,N
            KTP(J,I)=-1
            KTP(I,J)=KTPAIR(ANGL,P(1,I),P(1,J),KTP(J,I))
 200     CONTINUE
 210  CONTINUE
C---CALCULATE ALL SINGLE KT's
      DO 230 I=1,N
         KTS(I)=KTSING(ANGL,TYPE,P(1,I))
 230  CONTINUE
      KTMAX=0
C---MAIN LOOP
 300  CONTINUE
C---FIND MINIMUM MEMBER OF KTP
      CALL KTPMIN(KTP,NMAX,N,IMIN,JMIN)
C---FIND MINIMUM MEMBER OF KTS
      CALL KTSMIN(KTS,NMAX,N,KMIN)
C---STORE Y VALUE OF TRANSITION FROM N TO N-1 JETS
      KTMIN=KTP(IMIN,JMIN)
      IF ((TYPE.GE.2.AND.TYPE.LE.4).AND.
     &     (RSQ*KTS(KMIN).LE.KTMIN.OR.N.EQ.1))
     &     KTMIN=RSQ*KTS(KMIN)
      KT(N)=KTMIN
      Y(N)=KT(N)*ETSQ
C---IF MONO.GT.1, SEQUENCE IS SUPPOSED TO BE MONOTONIC, IF NOT, WARN
      IF (KTMIN.LT.KTMAX.AND.MONO.GT.1) CALL KTWARN('KTCLUS',1,*999)
      IF (KTMIN.GE.KTMAX) KTMAX=KTMIN
C---IF LOWEST KT IS TO A BEAM, THROW IT AWAY AND MOVE LAST ENTRY UP
      IF (KTMIN.EQ.RSQ*KTS(KMIN)) THEN
         CALL KTMOVE(P,KTP,KTS,NMAX,N,KMIN,1)
C---UPDATE HISTORY AND CROSS-REFERENCES
         HIST(N)=KMIN
         INJET(N)=KMIN
         DO 400 I=N,NN
            IF (INJET(I).EQ.KMIN) THEN
               KTLAST(I)=KTMAX
               INJET(I)=0
            ELSEIF (INJET(I).EQ.N) THEN
               INJET(I)=KMIN
            ENDIF
 400     CONTINUE
C---OTHERWISE MERGE JETS IMIN AND JMIN AND MOVE LAST ENTRY UP
      ELSE
         CALL KTMERG(P,KTP,KTS,NMAX,IMIN,JMIN,N,TYPE,ANGL,MONO,RECO)
         CALL KTMOVE(P,KTP,KTS,NMAX,N,JMIN,1)
C---UPDATE HISTORY AND CROSS-REFERENCES
         HIST(N)=IMIN*NMAX+JMIN
         INJET(N)=IMIN
         DO 600 I=N,NN
            IF (INJET(I).EQ.JMIN) THEN
               INJET(I)=IMIN
            ELSEIF (INJET(I).EQ.N) THEN
               INJET(I)=JMIN
            ENDIF
 600     CONTINUE
      ENDIF
C---THATS ALL THERE IS TO IT
      N=N-1
      IF (N.GT.1 .OR. N.GT.0.AND.(TYPE.GE.2.AND.TYPE.LE.4)) GOTO 300
      IF (N.EQ.1) THEN
         KT(N)=1E20
         Y(N)=KT(N)*ETSQ
      ENDIF
      RETURN
 999  RETURN 1
      END
C-----------------------------------------------------------------------
      FUNCTION KTPAIR(ANGL,P,Q,ANGLE)
      IMPLICIT NONE
C---CALCULATE LOCAL KT OF PAIR, USING ANGULAR SCHEME:
C   1=>ANGULAR, 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
C   WHERE f(eta,phi)=2(COSH(eta)-COS(phi)) IS THE QCD EMISSION METRIC
C---IF ANGLE<0, IT IS SET TO THE ANGULAR PART OF THE LOCAL KT ON RETURN
C   IF ANGLE>0, IT IS USED INSTEAD OF THE ANGULAR PART OF THE LOCAL KT
      INTEGER ANGL
      DOUBLE PRECISION P(9),Q(9),KTPAIR,R,KTMDPI,ANGLE,ETA,PHI,ESQ
C---COMPONENTS OF MOMENTA ARE PX,PY,PZ,E,1/P,PT,ETA,PHI,PT**2
      R=ANGLE
      IF (ANGL.EQ.1) THEN
         IF (R.LE.0) R=2*(1-(P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3))*(P(5)*Q(5)))
         ESQ=MIN(P(4),Q(4))**2
      ELSEIF (ANGL.EQ.2.OR.ANGL.EQ.3) THEN
         IF (R.LE.0) THEN
            ETA=P(7)-Q(7)
            PHI=KTMDPI(P(8)-Q(8))
            IF (ANGL.EQ.2) THEN
               R=ETA**2+PHI**2
            ELSE
               R=2*(COSH(ETA)-COS(PHI))
            ENDIF
         ENDIF
         ESQ=MIN(P(9),Q(9))
      ELSE
         CALL KTWARN('KTPAIR',200,*999)
         STOP
      ENDIF
      KTPAIR=ESQ*R
      IF (ANGLE.LT.0) ANGLE=R
 999  END
C-----------------------------------------------------------------------
      FUNCTION KTSING(ANGL,TYPE,P)
      IMPLICIT NONE
C---CALCULATE KT OF PARTICLE, USING ANGULAR SCHEME:
C   1=>ANGULAR, 2=>DeltaR, 3=>f(DeltaEta,DeltaPhi)
C---TYPE=1 FOR E+E-, 2 FOR EP, 3 FOR PE, 4 FOR PP
C   FOR EP, PROTON DIRECTION IS DEFINED AS -Z
C   FOR PE, PROTON DIRECTION IS DEFINED AS +Z
      INTEGER ANGL,TYPE
      DOUBLE PRECISION P(9),KTSING,COSTH,R,SMALL
      DATA SMALL/1E-4/
      IF (ANGL.EQ.1) THEN
         COSTH=P(3)*P(5)
         IF (TYPE.EQ.2) THEN
            COSTH=-COSTH
         ELSEIF (TYPE.EQ.4) THEN
            COSTH=ABS(COSTH)
         ELSEIF (TYPE.NE.1.AND.TYPE.NE.3) THEN
            CALL KTWARN('KTSING',200,*999)
            STOP
         ENDIF
         R=2*(1-COSTH)
C---IF CLOSE TO BEAM, USE APPROX 2*(1-COS(THETA))=SIN**2(THETA)
         IF (R.LT.SMALL) R=(P(1)**2+P(2)**2)*P(5)**2
         KTSING=P(4)**2*R
      ELSEIF (ANGL.EQ.2.OR.ANGL.EQ.3) THEN
         KTSING=P(9)
      ELSE
         CALL KTWARN('KTSING',201,*999)
         STOP
      ENDIF
 999  END
C-----------------------------------------------------------------------
      SUBROUTINE KTPMIN(A,NMAX,N,IMIN,JMIN)
      IMPLICIT NONE
C---FIND THE MINIMUM MEMBER OF A(NMAX,NMAX) WITH IMIN < JMIN <= N
      INTEGER NMAX,N,IMIN,JMIN,KMIN,I,J,K
C---REMEMBER THAT A(X+(Y-1)*NMAX)=A(X,Y)
C   THESE LOOPING VARIABLES ARE J=Y-2, I=X+(Y-1)*NMAX
      DOUBLE PRECISION A(*),AMIN
      K=1+NMAX
      KMIN=K
      AMIN=A(KMIN)
      DO 110 J=0,N-2
         DO 100 I=K,K+J
            IF (A(I).LT.AMIN) THEN
               KMIN=I
               AMIN=A(KMIN)
            ENDIF
 100     CONTINUE
         K=K+NMAX
 110  CONTINUE
      JMIN=KMIN/NMAX+1
      IMIN=KMIN-(JMIN-1)*NMAX
      END
C-----------------------------------------------------------------------
      SUBROUTINE KTSMIN(A,NMAX,N,IMIN)
      IMPLICIT NONE
C---FIND THE MINIMUM MEMBER OF A
      INTEGER N,NMAX,IMIN,I
      DOUBLE PRECISION A(NMAX)
      IMIN=1
      DO 100 I=1,N
         IF (A(I).LT.A(IMIN)) IMIN=I
 100  CONTINUE
      END
C-----------------------------------------------------------------------
      SUBROUTINE KTCOPY(A,N,B,ONSHLL)
      IMPLICIT NONE
C---COPY FROM A TO B. 5TH=1/(3-MTM), 6TH=PT, 7TH=ETA, 8TH=PHI, 9TH=PT**2
C   IF ONSHLL IS .TRUE. PARTICLE ENTRIES ARE PUT ON-SHELL BY SETTING E=P
      INTEGER I,N
      DOUBLE PRECISION A(4,N)
      LOGICAL ONSHLL
      DOUBLE PRECISION B(9,N),ETAMAX,SINMIN,EPS
      DATA ETAMAX,SINMIN,EPS/10,0,1E-6/
C---SINMIN GETS CALCULATED ON FIRST CALL
      IF (SINMIN.EQ.0) SINMIN=1/COSH(ETAMAX)
      DO 100 I=1,N
         B(1,I)=A(1,I)
         B(2,I)=A(2,I)
         B(3,I)=A(3,I)
         B(4,I)=A(4,I)
         B(5,I)=SQRT(A(1,I)**2+A(2,I)**2+A(3,I)**2)
         IF (ONSHLL) B(4,I)=B(5,I)
         IF (B(5,I).EQ.0) B(5,I)=1E-10
         B(5,I)=1/B(5,I)
         B(9,I)=A(1,I)**2+A(2,I)**2
         B(6,I)=SQRT(B(9,I))
         B(7,I)=B(6,I)*B(5,I)
         IF (B(7,I).GT.SINMIN) THEN
            B(7,I)=A(4,I)**2-A(3,I)**2
            IF (B(7,I).LE.EPS*B(4,I)**2.OR.ONSHLL) B(7,I)=B(9,I)
            B(7,I)=0.5*LOG((B(4,I)+ABS(B(3,I)))**2/B(7,I))
         ELSE
            B(7,I)=ETAMAX+2
         ENDIF
         B(7,I)=SIGN(B(7,I),B(3,I))
         IF (A(1,I).EQ.0 .AND. A(2,I).EQ.0) THEN
            B(8,I)=0
         ELSE
            B(8,I)=ATAN2(A(2,I),A(1,I))
         ENDIF
 100  CONTINUE
      END
C-----------------------------------------------------------------------
      SUBROUTINE KTMERG(P,KTP,KTS,NMAX,I,J,N,TYPE,ANGL,MONO,RECO)
      IMPLICIT NONE
C---MERGE THE Jth PARTICLE IN P INTO THE Ith PARTICLE
C   J IS ASSUMED GREATER THAN I. P CONTAINS N PARTICLES BEFORE MERGING.
C---ALSO RECALCULATING THE CORRESPONDING KTP AND KTS VALUES IF MONO.GT.0
C   FROM THE RECOMBINED ANGULAR MEASURES IF MONO.GT.1
C---NOTE THAT IF MONO.LE.0, TYPE AND ANGL ARE NOT USED
      INTEGER ANGL,RECO,TYPE,I,J,K,N,NMAX,MONO
      DOUBLE PRECISION P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX),PT,PTT,
     &     KTMDPI,KTUP,PI,PJ,ANG,KTPAIR,KTSING,ETAMAX,EPS
      KTUP(I,J)=KTP(MAX(I,J),MIN(I,J))
      DATA ETAMAX,EPS/10,1E-6/
      IF (J.LE.I) CALL KTWARN('KTMERG',200,*999)
C---COMBINE ANGULAR MEASURES IF NECESSARY
      IF (MONO.GT.1) THEN
         DO 100 K=1,N
            IF (K.NE.I.AND.K.NE.J) THEN
               IF (RECO.EQ.1) THEN
                  PI=P(4,I)
                  PJ=P(4,J)
               ELSEIF (RECO.EQ.2) THEN
                  PI=P(6,I)
                  PJ=P(6,J)
               ELSEIF (RECO.EQ.3) THEN
                  PI=P(9,I)
                  PJ=P(9,J)
               ELSE
                  CALL KTWARN('KTMERG',201,*999)
                  STOP
               ENDIF
               IF (PI.EQ.0.AND.PJ.EQ.0) THEN
                  PI=1
                  PJ=1
               ENDIF
               KTP(MAX(I,K),MIN(I,K))=
     &              (PI*KTUP(I,K)+PJ*KTUP(J,K))/(PI+PJ)
            ENDIF
 100     CONTINUE
      ENDIF
      IF (RECO.EQ.1) THEN
C---VECTOR ADDITION
         P(1,I)=P(1,I)+P(1,J)
         P(2,I)=P(2,I)+P(2,J)
         P(3,I)=P(3,I)+P(3,J)
         P(4,I)=P(4,I)+P(4,J)
         P(5,I)=SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
         IF (P(5,I).EQ.0) THEN
            P(5,I)=1
         ELSE
            P(5,I)=1/P(5,I)
         ENDIF
      ELSEIF (RECO.EQ.2) THEN
C---PT WEIGHTED ETA-PHI ADDITION
         PT=P(6,I)+P(6,J)
         IF (PT.EQ.0) THEN
            PTT=1
         ELSE
            PTT=1/PT
         ENDIF
         P(7,I)=(P(6,I)*P(7,I)+P(6,J)*P(7,J))*PTT
         P(8,I)=KTMDPI(P(8,I)+P(6,J)*PTT*KTMDPI(P(8,J)-P(8,I)))
         P(6,I)=PT
         P(9,I)=PT**2
      ELSEIF (RECO.EQ.3) THEN
C---PT**2 WEIGHTED ETA-PHI ADDITION
         PT=P(9,I)+P(9,J)
         IF (PT.EQ.0) THEN
            PTT=1
         ELSE
            PTT=1/PT
         ENDIF
         P(7,I)=(P(9,I)*P(7,I)+P(9,J)*P(7,J))*PTT
         P(8,I)=KTMDPI(P(8,I)+P(9,J)*PTT*KTMDPI(P(8,J)-P(8,I)))
         P(6,I)=P(6,I)+P(6,J)
         P(9,I)=P(6,I)**2
      ELSE
         CALL KTWARN('KTMERG',202,*999)
         STOP
      ENDIF
C---IF MONO.GT.0 CALCULATE NEW KT MEASURES. IF MONO.GT.1 USE ANGULAR ONES.
      IF (MONO.LE.0) RETURN
C---CONVERTING BETWEEN 4-MTM AND PT,ETA,PHI IF NECESSARY
      IF (ANGL.NE.1.AND.RECO.EQ.1) THEN
         P(9,I)=P(1,I)**2+P(2,I)**2
         P(7,I)=P(4,I)**2-P(3,I)**2
         IF (P(7,I).LE.EPS*P(4,I)**2) P(7,I)=P(9,I)
         IF (P(7,I).GT.0) THEN
            P(7,I)=0.5*LOG((P(4,I)+ABS(P(3,I)))**2/P(7,I))
            IF (P(7,I).GT.ETAMAX) P(7,I)=ETAMAX+2
         ELSE
            P(7,I)=ETAMAX+2
         ENDIF
         P(7,I)=SIGN(P(7,I),P(3,I))
         IF (P(1,I).NE.0.AND.P(2,I).NE.0) THEN
            P(8,I)=ATAN2(P(2,I),P(1,I))
         ELSE
            P(8,I)=0
         ENDIF
      ELSEIF (ANGL.EQ.1.AND.RECO.NE.1) THEN
         P(1,I)=P(6,I)*COS(P(8,I))
         P(2,I)=P(6,I)*SIN(P(8,I))
         P(3,I)=P(6,I)*SINH(P(7,I))
         P(4,I)=P(6,I)*COSH(P(7,I))
         IF (P(4,I).NE.0) THEN
            P(5,I)=1/P(4,I)
         ELSE
            P(5,I)=1
         ENDIF
      ENDIF
      ANG=0
      DO 200 K=1,N
         IF (K.NE.I.AND.K.NE.J) THEN
            IF (MONO.GT.1) ANG=KTUP(I,K)
            KTP(MIN(I,K),MAX(I,K))=
     &           KTPAIR(ANGL,P(1,I),P(1,K),ANG)
         ENDIF
 200  CONTINUE
      KTS(I)=KTSING(ANGL,TYPE,P(1,I))
 999  END
C-----------------------------------------------------------------------
      SUBROUTINE KTMOVE(P,KTP,KTS,NMAX,N,J,IOPT)
      IMPLICIT NONE
C---MOVE THE Nth PARTICLE IN P TO THE Jth POSITION
C---ALSO MOVING KTP AND KTS IF IOPT.GT.0
      INTEGER I,J,N,NMAX,IOPT
      DOUBLE PRECISION P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX)
      DO 100 I=1,9
         P(I,J)=P(I,N)
 100  CONTINUE
      IF (IOPT.LE.0) RETURN
      DO 110 I=1,J-1
         KTP(I,J)=KTP(I,N)
         KTP(J,I)=KTP(N,I)
 110  CONTINUE
      DO 120 I=J+1,N-1
         KTP(J,I)=KTP(I,N)
         KTP(I,J)=KTP(N,I)
 120  CONTINUE
      KTS(J)=KTS(N)
      END
C-----------------------------------------------------------------------
      FUNCTION KTMDPI(PHI)
      IMPLICIT NONE
C---RETURNS PHI, MOVED ONTO THE RANGE [-PI,PI)
      DOUBLE PRECISION KTMDPI,PHI,PI,TWOPI,THRPI,EPS
      PARAMETER (PI=3.141592654,TWOPI=6.283185307,THRPI=9.424777961)
      PARAMETER (EPS=1E-15)
      KTMDPI=PHI
      IF (KTMDPI.LE.PI) THEN
        IF (KTMDPI.GT.-PI) THEN
          GOTO 100
        ELSEIF (KTMDPI.GT.-THRPI) THEN
          KTMDPI=KTMDPI+TWOPI
        ELSE
          KTMDPI=-MOD(PI-KTMDPI,TWOPI)+PI
        ENDIF
      ELSEIF (KTMDPI.LE.THRPI) THEN
        KTMDPI=KTMDPI-TWOPI
      ELSE
        KTMDPI=MOD(PI+KTMDPI,TWOPI)-PI
      ENDIF
 100  IF (ABS(KTMDPI).LT.EPS) KTMDPI=0
      END
C-----------------------------------------------------------------------
      SUBROUTINE KTWARN(SUBRTN,ICODE,*)
C     DEALS WITH ERRORS DURING EXECUTION
C     SUBRTN = NAME OF CALLING SUBROUTINE
C     ICODE  = ERROR CODE:    - 99 PRINT WARNING & CONTINUE
C                          100-199 PRINT WARNING & JUMP
C                          200-    PRINT WARNING & STOP DEAD
C-----------------------------------------------------------------------
      INTEGER ICODE
      CHARACTER*6 SUBRTN
      WRITE (6,10) SUBRTN,ICODE
   10 FORMAT(/' KTWARN CALLED FROM SUBPROGRAM ',A6,': CODE =',I4/)
      IF (ICODE.LT.100) RETURN
      IF (ICODE.LT.200) RETURN 1
      STOP
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c
c spence function
c
      block data splint
      implicit double precision (a-h,o-z)
      common/spint/a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2
      data a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2/
     1 -0.250000000000000D0,
     2 -0.111111111111111D0,
     3 -0.010000000000000D0,
     4 -0.017006802721088D0,
     5 -0.019444444444444D0,
     6 -0.020661157024793D0,
     7 -0.021417300648069D0,
     8 -0.021948866377231D0,
     9 -0.022349233811171D0,
     1 -0.022663689135191D0,
     2  1.644934066848226D0/
      end
c
c spence function taking only real arguments
c
      function rsp(x)
      implicit double precision(a-h,o-z)
      common/spint/a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,zeta2
      x2=x*x
      if(x.gt.1.D0)then
        write(*,*)' argument greater than 1 passed to spence function'
        rsp=0.D0
        return
      endif
      if(x2.gt.1.D0.and.x.gt.0.5D0)then
        y=(x-1.D0)/x
        z=-log(1.D0-y)
        z2=z*z
        rsp=z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 +zeta2-log(x)*log(1.D0-x)+0.5D0*log(x)**2
        return
      elseif(x2.gt.1.D0.and.x.le.0.5D0)then
        y=1.D0/x
        z=-log(1.D0-y)
        z2=z*z
        rsp=-z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 -zeta2-0.5D0*log(-x)**2
        return
      elseif(x2.eq.1.D0)then
        rsp=zeta2
        return
      elseif(x2.le.1.D0.and.x.gt.0.5D0)then
        y=1.D0-x
        z=-log(1.D0-y)
        z2=z*z
        rsp=-z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
     3 +zeta2+z*log(1.D0-x)
       return
      elseif(x2.le.1.D0.and.x.le.0.5D0)then
        y=x
        z=-log(1.D0-y)
        z2=z*z
        rsp=z*(1.D0+a1*z*(1.D0+a2*z*(1.D0+a3*z2*(1.D0+a4*z2*
     1 (1.D0+a5*z2*(1.D0+a6*z2*(1.D0+a7*z2*(1.D0+a8*z2*(1.D0+a9*z2*
     2 (1.D0+a10*z2))))))))))
        return
      else
        write(*,*)' illegal x value in spence function'
        rsp=0.D0
      endif
      return
      end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE GFLIN1(ID,X,W)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NSIZE=200000,NMAX=2000)
      COMMON /GBOOK/ A(NSIZE)
      IF (ID.GT.NMAX) RETURN
      IS=INT(A(ID+2)+0.5)
      A(IS+9)=A(IS+9)+1.
      IOX=2
      IF(X.LT.A(IS+2)) IOX=1
      IF(X.GE.A(IS+3)) IOX=3
      A(IS+12+IOX)=A(IS+12+IOX)+W
      IF(IOX.NE.2) RETURN
      IX=INT((X-A(IS+2))/A(IS+4)+0.5)
      DX=(X-A(IS+2))/A(IS+4)+0.5-IX
      IF (IX.EQ.0) THEN
        A(IS+19+IX)=A(IS+19+IX)+W
      ELSEIF (IX.EQ.A(IS+1)) THEN
        A(IS+18+IX)=A(IS+18+IX)+W
      ELSE
        A(IS+18+IX)=A(IS+18+IX)+(1-DX)*W
        A(IS+19+IX)=A(IS+19+IX)+DX*W
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
