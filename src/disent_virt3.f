C-----------------------------------------------------------------------
C     Parity-violating part of the finite one-loop three-parton matrix
C     element, for Z and W exchange in DISENT's VIRTHR.
C
C     The one-loop amplitudes for 0 -> qbar q l lbar g are those of
C     Z. Bern, L. Dixon, D. Kosower, Nucl. Phys. B513 (1998) 3, as
C     implemented in MCFM by R.K. Ellis and J. Campbell (src/W1jet/A51.f,
C     A52.f, A5NLO.f, virt5.f, src/Need/spinoru.f, lnrat.f and
C     src/Wbb/lfunctions.f), from which this is a port (MCFM is GPL).
C     Momenta are all outgoing (negative energy for incoming particles)
C     in the layout (px,py,pz,E), as in DISENT.
C
C     The relation to DISENT (derivations/o2/oneloop/) was established
C     numerically: with all four helicity combinations summed, these
C     amplitudes reproduce DISENT's photon non-factorising one-loop part
C     (2*LEIV - EMSQ/2*ERTV)/EMSQ = CF*(CF*NX + CA*NY) exactly, via
C       I51 = C*(NX + 2*NY) + K51*T,   I52 = -C*NX + K52*T,
C     with C = -1/(8 pi^2), T the tree, I51/I52 the leading/subleading
C     colour interferences, and
C       K51 = -7/2 + pi^2/2 - (L13**2 + L23**2)/2,  K52 = 7/2 + L12**2/2,
C     L_ij = log(2|p_i.p_j|/EMSQ) (1 = incoming quark, 2 = outgoing quark,
C     3 = gluon). The parity-violating structures NX3, NY3 follow from the
C     same relations applied to the helicity combination
C     (LL + RR) - (LR + RL) (same minus opposite lepton-quark helicity).
C     The T-odd parts (odd under reflection of the event, sin(phi)
C     correlations from the absorptive parts) are dropped: they multiply
C     other coupling combinations and integrate to zero for
C     reflection-symmetric observables over DISENT's phase space.
C-----------------------------------------------------------------------
      SUBROUTINE VIRT3PV(P,NX3,NY3)
      IMPLICIT NONE
      DOUBLE PRECISION P(4,7),NX3,NY3,M5(4,5),S(5,5),MUSQ,DOT,PI,
     $     T(4),V51(4),V52(4),T3,I513,I523,L12,L13,L23,C,K51,K52
      DOUBLE COMPLEX ZA(5,5),ZB(5,5)
      INTEGER I,J,IP(5,4)
C---  helicities: LL, RR, LR, RL
      DATA IP/1,2,3,4,5, 2,1,4,3,5, 1,2,4,3,5, 2,1,3,4,5/
      PI=ATAN(1D0)*4
      MUSQ=-DOT(P,5,5)
      DO I=1,4
C---  qbar = incoming quark, q = outgoing quark, l = outgoing lepton,
C     lbar = incoming lepton, g = gluon
        M5(I,1)=-P(I,1)
        M5(I,2)= P(I,2)
        M5(I,3)= P(I,7)
        M5(I,4)=-P(I,6)
        M5(I,5)= P(I,3)
      ENDDO
      CALL V3SPIN(M5,ZA,ZB,S)
      DO J=1,4
        CALL V3HEL(IP(1,J),ZA,ZB,S,MUSQ,T(J),V51(J),V52(J))
      ENDDO
      T3  =T(1)+T(2)-T(3)-T(4)
      I513=V51(1)+V51(2)-V51(3)-V51(4)
      I523=V52(1)+V52(2)-V52(3)-V52(4)
      L12=LOG(2*ABS(DOT(P,1,2))/MUSQ)
      L13=LOG(2*ABS(DOT(P,1,3))/MUSQ)
      L23=LOG(2*ABS(DOT(P,2,3))/MUSQ)
      C=-1/(8*PI**2)
      K51=-3.5D0+PI**2/2-(L13**2+L23**2)/2
      K52= 3.5D0+L12**2/2
      NX3=-(I523-K52*T3)/C
      NY3=((I513-K51*T3)/C-NX3)/2
      END
C-----------------------------------------------------------------------
      SUBROUTINE V3HEL(IP,ZA,ZB,S,MUSQ,T,V51,V52)
      IMPLICIT NONE
C---  tree |A|^2 and Re(A_tree^* A_1loop) (leading colour V51, coefficient
C     of 1/N^2 V52), summed over the gluon helicities, for
C     0 -> qbar_R(IP1) q_L(IP2) l_L(IP3) lbar_R(IP4) g(IP5)  (MCFM virt5)
      INTEGER IP(5)
      DOUBLE PRECISION S(5,5),MUSQ,T,V51,V52
      DOUBLE COMPLEX ZA(5,5),ZB(5,5),LOM,A51M,A52M,LOP,A51P,A52P
      CALL V3A5(IP(1),IP(2),IP(3),IP(4),IP(5),ZA,ZB,S,MUSQ,LOM,A51M,
     $     A52M)
      CALL V3A5(IP(2),IP(1),IP(4),IP(3),IP(5),ZB,ZA,S,MUSQ,LOP,A51P,
     $     A52P)
      T=ABS(LOM)**2+ABS(LOP)**2
      V51=DBLE(DCONJG(LOM)*A51M)+DBLE(DCONJG(LOP)*A51P)
      V52=DBLE(DCONJG(LOM)*A52M)+DBLE(DCONJG(LOP)*A52P)
      END
C-----------------------------------------------------------------------
      SUBROUTINE V3A5(J1,J2,J3,J4,J5,ZA,ZB,S,MUSQ,LO,A51,A52)
      IMPLICIT NONE
C---  MCFM A5NLO: tree and one-loop (leading, subleading colour)
      INTEGER J1,J2,J3,J4,J5
      DOUBLE PRECISION S(5,5),MUSQ
      DOUBLE COMPLEX ZA(5,5),ZB(5,5),LO,A51,A52,V3A51,V3A52
      LO=-ZB(J1,J4)**2/(ZB(J2,J5)*ZB(J5,J1)*ZB(J4,J3))
      A51=V3A51(J2,J5,J1,J4,J3,ZB,ZA,S,MUSQ)
      A52=V3A52(J2,J1,J5,J4,J3,ZB,ZA,S,MUSQ)
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3A51(J1,J2,J3,J4,J5,ZA,ZB,S,MUSQ)
      IMPLICIT NONE
C---  MCFM A51 (leading colour), poles dropped
      INTEGER J1,J2,J3,J4,J5
      DOUBLE PRECISION S(5,5),MUSQ
      DOUBLE COMPLEX ZA(5,5),ZB(5,5),VCC,FCC,VSC,FSC,L12,L23,A5LOM,
     $     V3LNR,V3L0,V3L1,V3LSM1
      A5LOM=-ZA(J3,J4)**2/(ZA(J1,J2)*ZA(J2,J3)*ZA(J4,J5))
      L12=V3LNR(MUSQ,-S(J1,J2))
      L23=V3LNR(MUSQ,-S(J2,J3))
      VCC=-0.5D0*L12**2-0.5D0*L23**2-2*L23-4
      FCC=ZA(J3,J4)**2/(ZA(J1,J2)*ZA(J2,J3)*ZA(J4,J5))
     $     *(V3LSM1(-S(J1,J2),-S(J4,J5),-S(J2,J3),-S(J4,J5))
     $     -2*ZA(J3,J1)*ZB(J1,J5)*ZA(J5,J4)/ZA(J3,J4)
     $     *V3L0(-S(J2,J3),-S(J4,J5))/S(J4,J5))
      VSC=0.5D0*L23+1
      FSC=ZA(J3,J4)*ZA(J3,J1)*ZB(J1,J5)*ZA(J5,J4)
     $     /(ZA(J1,J2)*ZA(J2,J3)*ZA(J4,J5))*V3L0(-S(J2,J3),-S(J4,J5))
     $     /S(J4,J5)
     $     +0.5D0*(ZA(J3,J1)*ZB(J1,J5))**2*ZA(J4,J5)
     $     /(ZA(J1,J2)*ZA(J2,J3))*V3L1(-S(J2,J3),-S(J4,J5))/S(J4,J5)**2
      V3A51=(VCC+VSC)*A5LOM+FCC+FSC
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3A52(J1,J2,J3,J4,J5,ZA,ZB,S,MUSQ)
      IMPLICIT NONE
C---  MCFM A52 (subleading colour), poles dropped
      INTEGER J1,J2,J3,J4,J5
      DOUBLE PRECISION S(5,5),MUSQ
      DOUBLE COMPLEX ZA(5,5),ZB(5,5),VCC,FCC,VSC,FSC,L12,L45,A5LOM,
     $     V3LNR,V3L0,V3L1,V3LSM1
      L12=V3LNR(MUSQ,-S(J1,J2))
      L45=V3LNR(MUSQ,-S(J4,J5))
      A5LOM=ZA(J2,J4)**2/(ZA(J2,J3)*ZA(J3,J1)*ZA(J4,J5))
      VCC=-0.5D0*L12**2-2*L45-4
      FCC=-ZA(J2,J4)**2/(ZA(J2,J3)*ZA(J3,J1)*ZA(J4,J5))
     $     *V3LSM1(-S(J1,J2),-S(J4,J5),-S(J1,J3),-S(J4,J5))
     $     +ZA(J2,J4)*(ZA(J1,J2)*ZA(J3,J4)-ZA(J1,J4)*ZA(J2,J3))
     $     /(ZA(J2,J3)*ZA(J1,J3)**2*ZA(J4,J5))
     $     *V3LSM1(-S(J1,J2),-S(J4,J5),-S(J2,J3),-S(J4,J5))
     $     +2*ZB(J1,J3)*ZA(J1,J4)*ZA(J2,J4)/(ZA(J1,J3)*ZA(J4,J5))
     $     *V3L0(-S(J2,J3),-S(J4,J5))/S(J4,J5)
      VSC=0.5D0*L45+0.5D0
      FSC=ZA(J1,J4)**2*ZA(J2,J3)/(ZA(J1,J3)**3*ZA(J4,J5))
     $     *V3LSM1(-S(J1,J2),-S(J4,J5),-S(J2,J3),-S(J4,J5))
     $     -0.5D0*(ZA(J4,J1)*ZB(J1,J3))**2*ZA(J2,J3)
     $     /(ZA(J1,J3)*ZA(J4,J5))*V3L1(-S(J4,J5),-S(J2,J3))/S(J2,J3)**2
     $     +ZA(J1,J4)**2*ZA(J2,J3)*ZB(J3,J1)/(ZA(J1,J3)**2*ZA(J4,J5))
     $     *V3L0(-S(J4,J5),-S(J2,J3))/S(J2,J3)
     $     -ZA(J2,J1)*ZB(J1,J3)*ZA(J4,J3)*ZB(J3,J5)/ZA(J1,J3)
     $     *V3L1(-S(J4,J5),-S(J1,J2))/S(J1,J2)**2
     $     -ZA(J2,J1)*ZB(J1,J3)*ZA(J3,J4)*ZA(J1,J4)
     $     /(ZA(J1,J3)**2*ZA(J4,J5))*V3L0(-S(J4,J5),-S(J1,J2))/S(J1,J2)
     $     -0.5D0*ZB(J3,J5)*(ZB(J1,J3)*ZB(J2,J5)+ZB(J2,J3)*ZB(J1,J5))
     $     /(ZB(J1,J2)*ZB(J2,J3)*ZA(J1,J3)*ZB(J4,J5))
      V3A52=(VCC+VSC)*A5LOM+FCC+FSC
      END
C-----------------------------------------------------------------------
      SUBROUTINE V3SPIN(P,ZA,ZB,S)
      IMPLICIT NONE
C---  spinor products (MCFM spinoru), all momenta outgoing,
C     za(i,j)*zb(j,i) = s(i,j)
      DOUBLE PRECISION P(4,5),S(5,5),RT(5)
      DOUBLE COMPLEX ZA(5,5),ZB(5,5),C23(5),F(5)
      INTEGER I,J
      DO J=1,5
        ZA(J,J)=0
        ZB(J,J)=0
        S(J,J)=0
        IF (P(4,J).GT.0) THEN
          RT(J)=SQRT(P(4,J)+P(1,J))
          C23(J)=DCMPLX(P(3,J),-P(2,J))
          F(J)=1
        ELSE
          RT(J)=SQRT(-P(4,J)-P(1,J))
          C23(J)=DCMPLX(-P(3,J),P(2,J))
          F(J)=DCMPLX(0D0,1D0)
        ENDIF
      ENDDO
      DO I=2,5
        DO J=1,I-1
          S(I,J)=2*(P(4,I)*P(4,J)-P(1,I)*P(1,J)-P(2,I)*P(2,J)
     $         -P(3,I)*P(3,J))
          ZA(I,J)=F(I)*F(J)*(C23(I)*RT(J)/RT(I)-C23(J)*RT(I)/RT(J))
          ZB(I,J)=-S(I,J)/ZA(I,J)
          ZA(J,I)=-ZA(I,J)
          ZB(J,I)=-ZB(I,J)
          S(J,I)=S(I,J)
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3LNR(X,Y)
      IMPLICIT NONE
C---  log(x - i0) - log(y - i0)
      DOUBLE PRECISION X,Y,PI,TX,TY
      PI=ATAN(1D0)*4
      TX=0
      TY=0
      IF (X.LT.0) TX=1
      IF (Y.LT.0) TY=1
      V3LNR=DCMPLX(LOG(ABS(X/Y)),-PI*(TX-TY))
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3L0(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,D
      DOUBLE COMPLEX V3LNR
      D=1-X/Y
      IF (ABS(D).LT.1D-7) THEN
        V3L0=-1-D*(0.5D0+D/3)
      ELSE
        V3L0=V3LNR(X,Y)/D
      ENDIF
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3L1(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,D
      DOUBLE COMPLEX V3L0
      D=1-X/Y
      IF (ABS(D).LT.1D-7) THEN
        V3L1=-0.5D0-D/3*(1+0.75D0*D)
      ELSE
        V3L1=(V3L0(X,Y)+1)/D
      ENDIF
      END
C-----------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION V3LSM1(X1,Y1,X2,Y2)
      IMPLICIT NONE
C---  BDK Ls_{-1}; continuation checked against QCDLoop in all sign
C     regions (derivations/o2/oneloop)
      DOUBLE PRECISION X1,Y1,X2,Y2,R1,R2,O1,O2,DILOG,PISQO6
      DOUBLE COMPLEX D1,D2,V3LNR
      PISQO6=(ATAN(1D0)*4)**2/6
      R1=X1/Y1
      R2=X2/Y2
      O1=1-R1
      O2=1-R2
      IF (O1.GT.1) THEN
        D1=PISQO6-DILOG(R1)-V3LNR(X1,Y1)*LOG(O1)
      ELSE
        D1=DILOG(O1)
      ENDIF
      IF (O2.GT.1) THEN
        D2=PISQO6-DILOG(R2)-V3LNR(X2,Y2)*LOG(O2)
      ELSE
        D2=DILOG(O2)
      ENDIF
      V3LSM1=D1+D2+V3LNR(X1,Y1)*V3LNR(X2,Y2)-PISQO6
      END
