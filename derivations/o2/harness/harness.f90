! Evaluate DISENT's photon-exchange O(as^2) ingredients (src/libdisent.f)
! at phase-space points read from a file (kin.py layout), printing the
! individual structures of MATFOR (before external factors), CONTHR,
! ERTV and LEIV, for the point-by-point comparison with the FORM results.
! usage: harness <mode> <points file> [CF CA TR]
!   mode = matfor | conthr | virt
program harness
  implicit none
  integer :: SCHEME, NF
  double precision :: CF, CA, TR, PI, PISQ, HF, CUTOFF, EQ(-6:6), SCALE
  common /COLFAC/ CF, CA, TR, PI, PISQ, HF, CUTOFF, EQ, SCALE, SCHEME, NF
  double precision :: P(4,7), vin(28), EMSQ, A, B, C, D1, D2, E, QG, G, DOT
  double precision :: LEIA, LEIB, LEIC, LEID, LEIE, ERTA, ERTB, ERTC, ERTD, ERTE
  double precision :: NXS, NYS, NX3, NY3
  double precision :: ERTV, LEIV, V(4), VV, SC, cth, phi, C3PV, C3SY, X1PV, X2PV, X1SY, X2SY, DPV, DSY, EPV, ESY
  character(len=16) :: mode, arg
  integer :: IMODE
  common /HCOUP/ IMODE
  double precision :: M4(-6:6), SS4(-6:6), QM(4,7), S4(-6:6), JAC4, SSV
  integer :: J4, ifl, nok, iperm, perm3(3,6)
  double precision :: MS4(-6:6), PP(4,7)
  data perm3 /2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2/
  integer :: iunres, upos, kperm(6), lperm(6)
  character(len=8) :: envv
  ! K and L of SUBFOR's IPERM table
  data kperm /4,3,1,1,2,1/, lperm /1,1,4,3,1,2/
  character(len=256) :: fname
  integer :: i, ios
  external DOT, LEIA, LEIB, LEIC, LEID, LEIE, ERTA, ERTB, ERTC, ERTD, ERTE, ERTV, LEIV
  call get_command_argument(1, mode)
  call get_command_argument(2, fname)
  CF = 4d0/3; CA = 3; TR = 0.5d0
  if (command_argument_count() >= 5) then
     call get_command_argument(3, arg); read(arg,*) CF
     call get_command_argument(4, arg); read(arg,*) CA
     call get_command_argument(5, arg); read(arg,*) TR
  endif
  iunres = 0
  call get_environment_variable('UNRES', envv)
  if (len_trim(envv) > 0) read(envv,*) iunres
  IMODE = 0
  if (trim(mode) == 'subnc') IMODE = 1
  if (trim(mode) == 'subcc') IMODE = 2
  NF = 5; PI = atan(1d0)*4; PISQ = PI**2; HF = 0.5d0; CUTOFF = 1d-8
  EQ(0) = 0; EQ(1) = -1d0/3; EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo
  open(10, file=fname, status='old')
  do
     read(10, *, iostat=ios) vin
     if (ios /= 0) exit
     do i = 1, 7
        P(:,i) = vin(4*i-3:4*i)
     enddo
     EMSQ = -DOT(P,5,5)
     select case (trim(mode))
     case ('matfor')
        ! quark-initiated q -> q g g (as in MATFOR)
        A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4) &
             +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)+ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
        B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4) &
             +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)+ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))
        C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4) &
             +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)+ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
        QG = HF*(CF*A+(CF-CA/2)*B+CA*C)
        ! gluon-initiated g -> q qbar g
        A=2*(LEIA(P,P(1,6),2,3,-1,4)+LEIA(P,P(1,6),3,2,-1,4) &
             +LEIA(P,P(1,6),2,3,4,-1)+LEIA(P,P(1,6),3,2,4,-1))-EMSQ/2*( &
             ERTA(P,2,3,-1,4)+ERTA(P,3,2,-1,4)+ERTA(P,2,3,4,-1)+ERTA(P,3,2,4,-1))
        B=2*(LEIB(P,P(1,6),2,3,-1,4)+LEIB(P,P(1,6),3,2,-1,4) &
             +LEIB(P,P(1,6),2,3,4,-1)+LEIB(P,P(1,6),3,2,4,-1))-EMSQ/2*( &
             ERTB(P,2,3,-1,4)+ERTB(P,3,2,-1,4)+ERTB(P,2,3,4,-1)+ERTB(P,3,2,4,-1))
        C=2*(LEIC(P,P(1,6),2,3,-1,4)+LEIC(P,P(1,6),3,2,-1,4) &
             +LEIC(P,P(1,6),2,3,4,-1)+LEIC(P,P(1,6),3,2,4,-1))-EMSQ/2*( &
             ERTC(P,2,3,-1,4)+ERTC(P,3,2,-1,4)+ERTC(P,2,3,4,-1)+ERTC(P,3,2,4,-1))
        G = -(CF*A+(CF-CA/2)*B+CA*C)
        D1=2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*( &
             ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))
        D2=2*(LEID(P,P(1,6),-1,4,2,3)+LEID(P,P(1,6),2,3,-1,4))-EMSQ/2*( &
             ERTD(P,-1,4,2,3)+ERTD(P,2,3,-1,4))
        E=2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1) &
             +LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1) &
             +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4) &
             +LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4))-EMSQ/2*( &
             ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)+ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1) &
             +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)+ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
        write(*,'(6es25.16)') QG, G, D1, D2, E, EMSQ
     case ('conthr')
        ! V = a p2 + b p3 + ... perpendicular to the gluon (p3): use the
        ! vector of GTPERP-like construction from the point itself
        call vperp(P, 3, V)
        VV = V(4)**2-V(3)**2-V(2)**2-V(1)**2
        call CONTHR(P, V, VV, 2, -1, 3, SC)
        write(*,'(6es25.16)', advance='no') SC, VV, V
        call vperp(P, 1, V)
        VV = V(4)**2-V(3)**2-V(2)**2-V(1)**2
        call CONTHR(P, V, VV, 2, 3, -1, SC)
        write(*,'(6es25.16)') SC, VV, V
     case ('fo2')
        ! quark-initiated structures from MATFOR vs generated code (sym parts)
        A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4) &
             +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)+ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
        B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4) &
             +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)+ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))
        C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4) &
             +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)+ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
        QG = HF*(CF*A+(CF-CA/2)*B+CA*C)
        D1=2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*( &
             ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))
        E=2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1) &
             +LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1) &
             +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4) &
             +LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4))-EMSQ/2*( &
             ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)+ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1) &
             +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)+ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
        call FQGG3(P(1,6),P(1,1),P(1,2),P(1,3),P(1,4),X1PV,X2PV,X1SY,X2SY)
        call FD13(P(1,6),P(1,1),P(1,2),P(1,3),P(1,4),DPV,DSY)
        call FE3(P(1,6),P(1,1),P(1,2),P(1,3),P(1,4),EPV,ESY)
        write(*,'(3f20.15,3es14.5)') QG/((CF*X1SY+(CF-CA/2)*X2SY)/64), D1/(DSY/32), E/(-ESY/16), &
             (CF*X1PV+(CF-CA/2)*X2PV)/(CF*X1SY+(CF-CA/2)*X2SY), DPV/DSY, EPV/ESY
     case ('subph', 'subnc', 'subcc')
        ! symmetrise over the 6 labellings of the outgoing partons 2,3,4
        MS4 = 0; SS4 = 0; SSV = 1d6; nok = 0
        do iperm = 1, 6
           PP = P
           PP(:,2) = P(:,perm3(1,iperm)); PP(:,3) = P(:,perm3(2,iperm)); PP(:,4) = P(:,perm3(3,iperm))
           call MATFOR(PP, M4)
           MS4 = MS4 + M4
           ! label of the unresolved parton in this labelling
           upos = 0
           do J4 = 1, 3
              if (perm3(J4,iperm) == iunres) upos = J4 + 1
           enddo
           do J4 = 1, 6
              ! skip dipoles whose Born contains the unresolved parton (K or L)
              if (upos /= 0 .and. (kperm(J4) == upos .or. lperm(J4) == upos)) cycle
              call SUBFOR(J4, SSV, PP, QM, S4, JAC4, *200)
              SS4 = SS4 + S4; nok = nok + 1
200           continue
           enddo
        enddo
        M4 = MS4
        write(*,'(i3,6es16.7)') nok, (M4(ifl)/SS4(ifl), ifl=-2,2), -DOT(P,5,5)
     case ('soft')
        ! MATFOR quark-channel q->qgg part vs eikonal(gluon 4) x MATTHR Born
        A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4) &
             +LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)+ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
        B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4) &
             +LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)+ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))
        C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4) &
             +LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3))-EMSQ/2*( &
             ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)+ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
        QG = HF*(CF*A+(CF-CA/2)*B+CA*C)*256*PI**4*CF/EMSQ*(4*PI/137)**2*4/EMSQ
        PP = P; PP(:,4) = 0
        call MATTHR(PP, M4)
        ! eikonal: -16 pi^2 sum_{i<j} Ti.Tj (pi.pj)/((pi.p4)(pj.p4)), crossed p1 -> -p1
        ! T1.T2 = CA/2-CF, T1.T3 = T2.T3 = -CA/2 (1 = incoming quark as outgoing antiquark)
        SC = -16*PISQ*( (CA/2-CF)*(-DOT(P,1,2))/((-DOT(P,1,4))*DOT(P,2,4)) &
             - CA/2*(-DOT(P,1,3))/((-DOT(P,1,4))*DOT(P,3,4)) &
             - CA/2*DOT(P,2,3)/(DOT(P,2,4)*DOT(P,3,4)) )
        write(*,'(3es20.10)') QG*EQ(2)**2/(HF*SC*M4(2)), QG, M4(2)
     case ('softdip')
        ! per-dipole SUBFOR (quark I=2) vs eikonal pieces x Born, soft gluon 4
        PP = P; PP(:,4) = 0
        call MATTHR(PP, M4)
        write(*,'(a,3es14.5)') ' eik x HF x Born: (12) (13) (23):', &
             -16*PISQ*(CA/2-CF)*(-DOT(P,1,2))/((-DOT(P,1,4))*DOT(P,2,4))*HF*M4(2), &
             -16*PISQ*(-CA/2)*(-DOT(P,1,3))/((-DOT(P,1,4))*DOT(P,3,4))*HF*M4(2), &
             -16*PISQ*(-CA/2)*DOT(P,2,3)/(DOT(P,2,4)*DOT(P,3,4))*HF*M4(2)
        SSV = 1d6
        do J4 = 1, 6
           S4 = 0
           call SUBFOR(J4, SSV, P, QM, S4, JAC4, *201)
201        continue
           write(*,'(a,i2,es14.5)') '  perm', J4, S4(2)
        enddo
     case ('fcc')
        call FEXX3(P(1,6),P(1,1),P(1,2),P(1,3),P(1,4),EPV,ESY)
        write(*,'(2es25.16)',advance='no') ESY, EPV
        call FEXY3(P(1,6),P(1,1),P(1,2),P(1,3),P(1,4),EPV,ESY)
        write(*,'(2es25.16)') ESY, EPV
     case ('virt3')
        call virt3sy(P, NXS, NYS)
        call VIRT3PV(P, NX3, NY3)
        ! DISENT (2 LEIV - EMSQ/2 ERTV)/EMSQ for the current CF, CA
        write(*,'(4es25.16)') (2*LEIV(P,P(1,6),2,-1,3)-EMSQ/2*ERTV(P,2,-1,3))/EMSQ, &
             CF*(CF*NXS+CA*NYS), NX3, NY3
     case ('conthr3')
        call vperp(P, 3, V)
        VV = V(4)**2-V(3)**2-V(2)**2-V(1)**2
        call CONTHR(P, V, VV, 2, -1, 3, SC)
        call FCONTH3(P(1,6), P(1,1), P(1,2), P(1,3), V, C3PV, C3SY)
        write(*,'(4es25.16)') SC, C3SY/VV/EMSQ**2, C3PV/VV/EMSQ**2, SC/(C3SY/VV/EMSQ**2)
     case ('virtee')
        write(*,'(5es25.16)') ERTV(P,1,2,3), LEIV(P,P(1,6),1,2,3), &
             ERTV(P,1,2,3), LEIV(P,P(1,6),1,2,3), DOT(P,5,5)
     case ('virt')
        write(*,'(5es25.16)') ERTV(P,2,-1,3), LEIV(P,P(1,6),2,-1,3), &
             ERTV(P,2,3,-1), LEIV(P,P(1,6),2,3,-1), EMSQ
     end select
  enddo
contains
  ! a space-like vector perpendicular to P(,ig), built from P(,6)
  subroutine vperp(P, ig, V)
    double precision, intent(in) :: P(4,7)
    integer, intent(in) :: ig
    double precision, intent(out) :: V(4)
    double precision :: a, b
    a = DOT(P,6,ig); b = DOT(P,2,ig)
    V = P(:,6) - a/b * P(:,2)
  end subroutine vperp
end program harness

! symmetric (photon) version of VIRT3PV, for testing
subroutine virt3sy(P, NX, NY)
  implicit none
  double precision :: P(4,7), NX, NY, M5(4,5), S(5,5), MUSQ, DOT, PI, T(4), V51(4), V52(4)
  double precision :: TS, I51, I52, L12, L13, L23, C, K51, K52
  double complex :: ZA(5,5), ZB(5,5)
  integer :: I, J, IP(5,4)
  data IP /1,2,3,4,5, 2,1,4,3,5, 1,2,4,3,5, 2,1,3,4,5/
  external DOT
  PI = atan(1d0)*4; MUSQ = -DOT(P,5,5)
  do I = 1, 4
     M5(I,1) = -P(I,1); M5(I,2) = P(I,2); M5(I,3) = P(I,7); M5(I,4) = -P(I,6); M5(I,5) = P(I,3)
  enddo
  call V3SPIN(M5, ZA, ZB, S)
  do J = 1, 4
     call V3HEL(IP(1,J), ZA, ZB, S, MUSQ, T(J), V51(J), V52(J))
  enddo
  TS = sum(T); I51 = sum(V51); I52 = sum(V52)
  L12 = log(2*abs(DOT(P,1,2))/MUSQ); L13 = log(2*abs(DOT(P,1,3))/MUSQ); L23 = log(2*abs(DOT(P,2,3))/MUSQ)
  C = -1/(8*PI**2); K51 = -3.5d0 + PI**2/2 - (L13**2+L23**2)/2; K52 = 3.5d0 + L12**2/2
  NX = -(I52 - K52*TS)/C
  NY = ((I51 - K51*TS)/C - NX)/2
end subroutine virt3sy

! couplings for DISENT (called from libdisent.f): photon (IMODE 0), a
! generic parity-violating NC set (1: c2 = 1, c3 = +-0.6), or a pure
! left-handed CC set (2: u-type and dbar-type with c2 = 1, c3 = +-1)
subroutine hsetcoup(eq, c2n, c3n, c2c, c3c)
  implicit none
  double precision, intent(in) :: eq(-6:6)
  double precision, intent(out) :: c2n(-6:6), c3n(-6:6), c2c(-6:6), c3c(-6:6)
  integer :: IMODE, i
  common /HCOUP/ IMODE
  c2n = 0; c3n = 0; c2c = 0; c3c = 0
  if (IMODE == 0) then
     c2n = eq**2
  elseif (IMODE == 1) then
     do i = 1, 5
        c2n(i) = 1; c2n(-i) = 1; c3n(i) = 0.6d0; c3n(-i) = -0.6d0
     enddo
  else
     do i = 2, 4, 2
        c2c(i) = 1; c3c(i) = 1; c2c(-(i-1)) = 1; c3c(-(i-1)) = -1
     enddo
  endif
end subroutine hsetcoup

subroutine disent_couplings(Q2, eq, c2, c3, cg)
  implicit none
  double precision, intent(in)  :: Q2, eq(-6:6)
  double precision, intent(out) :: c2(-6:6), c3(-6:6), cg(6)
  double precision :: c2n(-6:6), c3n(-6:6), c2c(-6:6), c3c(-6:6)
  integer :: i
  call hsetcoup(eq, c2n, c3n, c2c, c3c)
  c2 = c2n + c2c; c3 = c3n + c3c
  do i = 1, 6
     cg(i) = (c2(i) + c2(-i))/2
  enddo
end subroutine disent_couplings

subroutine disent_couplings4(Q2, eq, c2n, c3n, c2c, c3c, cg)
  implicit none
  double precision, intent(in)  :: Q2, eq(-6:6)
  double precision, intent(out) :: c2n(-6:6), c3n(-6:6), c2c(-6:6), c3c(-6:6), cg(6)
  integer :: i
  call hsetcoup(eq, c2n, c3n, c2c, c3c)
  do i = 1, 6
     cg(i) = (c2n(i) + c2c(i) + c2n(-i) + c2c(-i))/2
  enddo
end subroutine disent_couplings4
