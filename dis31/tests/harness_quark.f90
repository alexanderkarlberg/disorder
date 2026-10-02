! Harness: incoming quark, photon exchange. MCFM crossed: gamma* q -> q g g
! (z2jetsq) + gamma* q -> q Q Qbar for all five flavours Q (ampqqb_qqb),
! including the identical-quark interference (both signs tried), against
! DISENT's MATFOR M(i); both summed over the labellings of the outgoing
! partons 2, 3, 4. The ratio must be constant (and equal for i = 1, 2 up
! to the known factors).
program harness31q
  use types, only: dp
  use mod_parameters
  implicit none
  integer :: SCHEME_D, NF
  double precision :: CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ(-6:6), SCALE_D
  common /COLFAC/ CF, CA, TR, PI_D, PISQ, HF, CUTOFF_D, EQ, SCALE_D, SCHEME_D, NF
  integer, parameter :: mxpart = 14
  complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
  common /zprods/ za, zb
  real(dp), parameter :: xn = 3
  real(dp) :: P(4,7), PP(4,7), M4(-6:6), r(8), pm(mxpart,4), msq(2,2)
  real(dp) :: dsum(2), msum(2,2), rat(2,2,20), cq, cpair, t4, tid(2)
  real(dp) :: dp4(4), dpc(4), mp(4), mpa(4), pr(4,20), ev(4), evs(4)
  integer :: kx, hx, hy
  complex(dp) :: Adir(2,2,2), Bdir(2,2,2), Aexc(2,2,2), Bexc(2,2,2)
  real(dp) :: emc, emcs
  integer, parameter :: swp(2) = [2, 1]
  complex(dp) :: Aq(2,2,2), Bq(2,2,2), Ax(2,2,2), Bx(2,2,2), aa, bb
  integer :: perm(3,6), ip, ipt, i, iq, ipair, j1, j2, j3, is
  perm = reshape([2,3,4, 2,4,3, 3,2,4, 3,4,2, 4,2,3, 4,3,2], [3,6])

  call set_parameters()
  CF = 4.0_dp/3.0_dp; CA = 3; TR = 0.5_dp; NF = 5
  PI_D = atan(1d0)*4; PISQ = PI_D**2; HF = 0.5_dp; CUTOFF_D = 1d-8
  SCHEME_D = 0; SCALE_D = 1
  EQ(0) = 0; EQ(1) = -1D0/3; EQ(2) = EQ(1) + 1
  do i = 1, 6
     if (i > 2) EQ(i) = EQ(i-2)
     EQ(-i) = -EQ(i)
  enddo

  do ipt = 1, 20
     call random_number(r)
     call dis_point(r, P)
     dsum = 0; msum = 0; dpc = 0; mp = 0; evs = 0; emcs = 0
     do ip = 1, 6
        PP = P
        PP(:,2) = P(:,perm(1,ip)); PP(:,3) = P(:,perm(2,ip)); PP(:,4) = P(:,perm(3,ip))
        call MATFOR(PP, M4)
        dsum = dsum + M4(1:2)
        call dpieces(PP, dp4)
        dpc = dpc + dp4
        ! slots: 1 = -q in, 2 = -electron in, 3 = electron out, 4 = q out,
        ! 5, 6 = the other two outgoing partons (g g, or Q Qbar)
        pm = 0
        ! slots (MCFM's four-quark routine has the leptons at 3, 4):
        ! 1 = -q in, 2 = q out, 3 = electron out, 4 = -electron in, 5, 6 = the
        ! other two outgoing partons (g g, or Q Qbar)
        pm(1,:) = -P(:,1); pm(4,:) = -P(:,6); pm(3,:) = P(:,7)
        pm(2,:) = P(:,perm(1,ip)); pm(5,:) = P(:,perm(2,ip)); pm(6,:) = P(:,perm(3,ip))
        call spinoru(6, pm, za, zb)
        call z2jetsq(2, 1, 3, 4, 5, 6, za, zb, msq)
        call ampqqb_qqb(2, 1, 5, 6, Aq, Bq)        ! q line (2,1), Q line (5,6)
        call ampqqb_qqb(5, 1, 2, 6, Ax, Bx)      ! identical: the two quarks 2, 5 exchanged
        call ampqqb_qqb(1, 2, 6, 5, Adir, Bdir)  ! MCFM's identical-quark construction (crossed)
        call ampqqb_qqb(1, 5, 2, 6, Aexc, Bexc)
        ! MCFM pieces, unit charges: q g g, |A|^2 (boson on the incoming line),
        ! |B|^2 (on the pair), identical-quark interference
        mpa = 0
        mpa(1) = (xn/4)*0.5_dp*sum(msq)
        do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
           mpa(2) = mpa(2) + 4*abs(Aq(j1,j2,j3))**2
           mpa(3) = mpa(3) + 4*abs(Bq(j1,j2,j3))**2
           if (j1 == j2) mpa(4) = mpa(4) + 0.5_dp*4*2/xn*real((Aq(j1,j2,j3) + Bq(j1,j2,j3))*conjg(Ax(j1,j2,j3) + Bx(j1,j2,j3)), dp)
        enddo; enddo; enddo
        mp = mp + mpa
        ! identical-quark interference, helicity pairings of the exchanged amplitude
        do kx = 1, 4
           ev(kx) = 0
           do j1 = 1, 2; do j3 = 1, 2
              hx = merge(j1, 3 - j1, kx <= 2); hy = merge(j1, 3 - j1, mod(kx,2) == 1)
              ev(kx) = ev(kx) + 0.5_dp*4*2/xn*real((Aq(j1,j1,j3) + Bq(j1,j1,j3))*conjg(Ax(hx,hy,j3) + Bx(hx,hy,j3)), dp)
           enddo; enddo
        enddo
        evs = evs + ev
        ! MCFM's identical q q -> q q construction, crossed: MCFM slots (1,2,5,6)
        ! = ours (1,6,2,5); direct ampqqb_qqb(1,5,2,6) with j2 swapped and B
        ! negated, exchange ampqqb_qqb(1,6,5,2)
        emc = 0
        do j1 = 1, 2; do j3 = 1, 2
           emc = emc + 0.5_dp*4*2/xn*real((Adir(j1,swp(j1),j3) - Bdir(j1,swp(j1),j3))*conjg(Aexc(j1,j1,j3) + Bexc(j1,j1,j3)), dp)
        enddo; enddo
        emcs = emcs + emc
        do iq = 1, 2
           cq = EQ(iq)
           ! q g g (identical gluons: 1/2)
           msum(iq,:) = msum(iq,:) + cq**2*(xn/4)*0.5_dp*sum(msq)
           do ipair = 1, 5
              cpair = EQ(ipair)
              if (ipair /= iq) then
                 t4 = 0
                 do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
                    t4 = t4 + abs(cq*Aq(j1,j2,j3) + cpair*Bq(j1,j2,j3))**2
                 enddo; enddo; enddo
                 msum(iq,:) = msum(iq,:) + 4*t4
              else
                 ! identical quarks, MCFM's construction (direct with j2 swapped
                 ! and B negated, exchange as called); photon charges equal
                 tid = 0
                 do j1 = 1, 2; do j2 = 1, 2; do j3 = 1, 2
                    aa = cq*(Adir(j1,swp(j2),j3) - Bdir(j1,swp(j2),j3))
                    bb = cq*(Aexc(j1,j2,j3) + Bexc(j1,j2,j3))
                    tid = tid + cq**2*(abs(Aq(j1,j2,j3))**2 + abs(Bq(j1,j2,j3))**2 + abs(Ax(j1,j2,j3))**2 + abs(Bx(j1,j2,j3))**2)
                    if (j2 == j1) tid = tid + 2/xn*real(aa*conjg(cq*(Aexc(j1,j1,j3) + Bexc(j1,j1,j3))), dp)
                 enddo; enddo; enddo
                 msum(iq,:) = msum(iq,:) + 0.5_dp*4*tid
              endif
           enddo
        enddo
     enddo
     do iq = 1, 2
        rat(iq,:,ipt) = dsum(iq)/msum(iq,:)
     enddo
     pr(:,ipt) = dpc/mp
     write(*,'(a,i3,a,2es20.12)') ' recon ', ipt, ' MATFOR / pieces (d, u):', &
          & dsum(1)/(EQ(1)**2*(dpc(1) + NF*TR*dpc(2) + dpc(4)) + sum(EQ(1:5)**2)*TR*dpc(3)), &
          & dsum(2)/(EQ(2)**2*(dpc(1) + NF*TR*dpc(2) + dpc(4)) + sum(EQ(1:5)**2)*TR*dpc(3))
     write(*,'(a,i3,a,4es16.8)') ' pieces', ipt, ' DISENT/MCFM (qgg, D1, D2, E):', pr(:,ipt)
     write(*,'(a,i3,a,es16.8)') ' EMCFM ', ipt, ' DISENT E / MCFM-construction interference:', dpc(4)/emcs
     write(*,'(a,i3,a,4es16.8)') ' Epair ', ipt, ' DISENT E / MCFM interference (hh, h-h, -hh, -h-h):', dpc(4)/evs
     write(*,'(a,i3,a,2es12.4,a,2es20.12)') ' point', ipt, '  DISENT d, u', dsum, &
          & '  ratio d (interference +, -)', rat(1,:,ipt)
  enddo
  do iq = 1, 2
     do is = 1, 2
        write(*,'(a,i2,a,i2,a,es20.12,a,es10.2)') ' flavour', iq, ' sign', 3 - 2*is, ': mean ratio', &
             & sum(rat(iq,is,:))/20, '  spread', (maxval(rat(iq,is,:)) - minval(rat(iq,is,:)))/abs(sum(rat(iq,is,:))/20)
     enddo
  enddo

contains

  ! DISENT's photon-exchange pieces of MATFOR (unnormalised): q g g, D1, D2, E
  subroutine dpieces(P, d)
    real(dp), intent(in) :: P(4,7)
    real(dp), intent(out) :: d(4)
    real(dp), external :: LEIA, LEIB, LEIC, LEID, LEIE, ERTA, ERTB, ERTC, ERTD, ERTE, DOT
    real(dp) :: A, B, C, EMSQ
    EMSQ = -DOT(P,5,5)
    A=2*(LEIA(P,P(1,6),-1,2,3,4)+LEIA(P,P(1,6),2,-1,3,4)+LEIA(P,P(1,6),-1,2,4,3)+LEIA(P,P(1,6),2,-1,4,3)) &
         & -EMSQ/2*(ERTA(P,-1,2,3,4)+ERTA(P,2,-1,3,4)+ERTA(P,-1,2,4,3)+ERTA(P,2,-1,4,3))
    B=2*(LEIB(P,P(1,6),-1,2,3,4)+LEIB(P,P(1,6),2,-1,3,4)+LEIB(P,P(1,6),-1,2,4,3)+LEIB(P,P(1,6),2,-1,4,3)) &
         & -EMSQ/2*(ERTB(P,-1,2,3,4)+ERTB(P,2,-1,3,4)+ERTB(P,-1,2,4,3)+ERTB(P,2,-1,4,3))
    C=2*(LEIC(P,P(1,6),-1,2,3,4)+LEIC(P,P(1,6),2,-1,3,4)+LEIC(P,P(1,6),-1,2,4,3)+LEIC(P,P(1,6),2,-1,4,3)) &
         & -EMSQ/2*(ERTC(P,-1,2,3,4)+ERTC(P,2,-1,3,4)+ERTC(P,-1,2,4,3)+ERTC(P,2,-1,4,3))
    d(1) = HF*(CF*A+(CF-CA/2)*B+CA*C)
    d(2) = 2*(LEID(P,P(1,6),4,-1,3,2)+LEID(P,P(1,6),3,2,4,-1))-EMSQ/2*(ERTD(P,4,-1,3,2)+ERTD(P,3,2,4,-1))
    d(3) = 2*(LEID(P,P(1,6),-1,4,2,3)+LEID(P,P(1,6),2,3,-1,4))-EMSQ/2*(ERTD(P,-1,4,2,3)+ERTD(P,2,3,-1,4))
    d(4) = 2*(LEIE(P,P(1,6),4,-1,3,2)+LEIE(P,P(1,6),3,2,4,-1)+LEIE(P,P(1,6),4,-1,2,3)+LEIE(P,P(1,6),2,3,4,-1) &
         & +LEIE(P,P(1,6),-1,4,2,3)+LEIE(P,P(1,6),2,3,-1,4)+LEIE(P,P(1,6),-1,4,3,2)+LEIE(P,P(1,6),3,2,-1,4)) &
         & -EMSQ/2*(ERTE(P,4,-1,3,2)+ERTE(P,3,2,4,-1)+ERTE(P,4,-1,2,3)+ERTE(P,2,3,4,-1) &
         & +ERTE(P,-1,4,2,3)+ERTE(P,2,3,-1,4)+ERTE(P,-1,4,3,2)+ERTE(P,3,2,-1,4))
    d(4) = HF*(CF-CA/2)*d(4)
    d = d/EMSQ**2   ! MATFOR's final normalisation is proportional to 1/EMSQ^2
  end subroutine dpieces

  subroutine dis_point(r, P)
    real(dp), intent(in) :: r(8)
    real(dp), intent(out) :: P(4,7)
    real(dp) :: Q, y, xi, E, E1, W, beta, m23, pa(4), pb(4), pc(4), X(4), d(3)
    Q = 10 + 90*r(1); y = 0.1_dp + 0.8_dp*r(2); xi = 0.05_dp + 0.8_dp*r(3)
    E = Q/2; E1 = E/xi; W = sqrt(4*E*E1 - 4*E**2); beta = (E1 - 2*E)/E1
    P = 0
    P(:,1) = [0.0_dp, 0.0_dp, E1, E1]
    P(:,5) = [0.0_dp, 0.0_dp, -2*E, 0.0_dp]
    P(:,6) = [E/y*2*sqrt(1-y), 0.0_dp, -E, E/y*(2-y)]
    P(:,7) = [E/y*2*sqrt(1-y), 0.0_dp,  E, E/y*(2-y)]
    m23 = W*(0.1_dp + 0.8_dp*r(4))
    d = unitv(2*r(5) - 1, 6.283185307_dp*r(6))
    pa = [(W*W - m23*m23)/(2*W)*d, (W*W - m23*m23)/(2*W)]
    X = [-pa(1:3), W - pa(4)]
    call twobody(X, unitv(2*r(7) - 1, 6.283185307_dp*r(8)), pb, pc)
    P(:,2) = boostz(pa, beta); P(:,3) = boostz(pb, beta); P(:,4) = boostz(pc, beta)
  end subroutine dis_point

  function unitv(cth, phi) result(d)
    real(dp), intent(in) :: cth, phi
    real(dp) :: d(3)
    d = [sqrt(1-cth**2)*cos(phi), sqrt(1-cth**2)*sin(phi), cth]
  end function unitv

  subroutine twobody(Ptot, dir, pa, pb)
    real(dp), intent(in)  :: Ptot(4), dir(3)
    real(dp), intent(out) :: pa(4), pb(4)
    real(dp) :: Mm, h
    Mm = sqrt(Ptot(4)**2 - sum(Ptot(1:3)**2)); h = Mm/2
    pa = boostv([h*dir, h], Ptot); pb = boostv([-h*dir, h], Ptot)
  end subroutine twobody

  function boostv(p, Pt) result(q)
    real(dp), intent(in) :: p(4), Pt(4)
    real(dp) :: q(4), Mm, bp, f
    Mm = sqrt(Pt(4)**2 - sum(Pt(1:3)**2))
    bp = dot_product(p(1:3), Pt(1:3))
    q(4) = (Pt(4)*p(4) + bp)/Mm
    f = (bp/(Pt(4) + Mm) + p(4))/Mm
    q(1:3) = p(1:3) + f*Pt(1:3)
  end function boostv

  function boostz(p, beta) result(q)
    real(dp), intent(in) :: p(4), beta
    real(dp) :: q(4), g
    g = 1/sqrt(1 - beta**2)
    q = p
    q(3) = g*(p(3) + beta*p(4)); q(4) = g*(p(4) + beta*p(3))
  end function boostz
end program harness31q
