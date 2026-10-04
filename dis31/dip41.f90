!-----------------------------------------------------------------------
! Catani-Seymour dipoles for DIS 4+1 -> 3+1 (photon exchange), one
! incoming parton (Catani, Seymour, hep-ph/9605323, section 5, epsilon = 0):
!   FF  D_{ij,k}: emitter pair i, j and spectator k in the final state,
!   FI  D_ij^a:   final-state pair, the incoming parton as spectator,
!   IF  D^{ai}_k: the incoming parton a and final-state i, spectator k.
! (II dipoles do not occur: the lepton is colourless.)
!
! Layouts as me41 (P4(4,8): 1 incoming parton, 2-5 outgoing, 6 q, 7, 8
! leptons) and born31 (P3(4,7)); normalisation as me41, i.e. each 8 pi
! alpha_s of CS becomes 16 pi^2 (matrix elements divided by
! (alpha_s/2pi)^n). Colour and spin correlations from born31; the Born
! averages are those of born31 (the averaged matrix elements of the mapped
! Born), so that sum_dipoles -> me41 in every single-unresolved limit.
!
! dip41_list(P4, fl4, nd, P3, fl3, val): all nd dipoles, each with its
!   mapped Born momenta P3(:,:,id), flavours fl3(:,id) and value val(id).
!-----------------------------------------------------------------------
module dip41
  use born31
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: CF = 4.0_dp/3, CA = 3, TR = 0.5_dp
  real(dp), parameter :: c16 = 16*pi**2
  integer, parameter, public :: dip41_max = 48
  ! labels of the last dip41_list call (diagnostics): type (1 FF, 2 FI, 3 IF),
  ! emitter i (IF: the final-state parton i), j (FF/FI) and spectator k
  integer, public, save :: dip41_lab(4,dip41_max)
  public :: dip41_list
contains

  subroutine dip41_list(P4, fl4, nd, P3, fl3, val)
    real(dp), intent(in) :: P4(4,8)
    integer, intent(in) :: fl4(5)
    integer, intent(out) :: nd
    real(dp), intent(out) :: P3(4,7,dip41_max), val(dip41_max)
    integer, intent(out) :: fl3(4,dip41_max)
    integer :: i, j, k, fij, fa
    nd = 0
    ! final-state pairs (unordered; for q g the quark is i)
    do i = 2, 5
       do j = 2, 5
          if (j == i) cycle
          if (.not. ff_pair(fl4(i), fl4(j), i, j, fij)) cycle
          do k = 1, 5
             if (k == i .or. k == j) cycle
             nd = nd + 1
             dip41_lab(:,nd) = [merge(2, 1, k == 1), i, j, k]
             if (k == 1) then
                call dip_fi(P4, fl4, i, j, fij, P3(:,:,nd), fl3(:,nd), val(nd))
             else
                call dip_ff(P4, fl4, i, j, k, fij, P3(:,:,nd), fl3(:,nd), val(nd))
             endif
          enddo
       enddo
    enddo
    ! initial-state emitter a = 1 with final-state i
    do i = 2, 5
       if (.not. if_pair(fl4(1), fl4(i), fa)) cycle
       do k = 2, 5
          if (k == i) cycle
          nd = nd + 1
          dip41_lab(:,nd) = [3, i, 0, k]
          call dip_if(P4, fl4, i, k, fa, P3(:,:,nd), fl3(:,nd), val(nd))
       enddo
    enddo
  end subroutine dip41_list

  ! can final-state partons i, j (flavours fi, fj) be merged? each unordered
  ! pair once: q g with the quark first, g g and q qbar with i < j
  logical function ff_pair(fi, fj, i, j, fij)
    integer, intent(in) :: fi, fj, i, j
    integer, intent(out) :: fij
    ff_pair = .false.; fij = 0
    if (fi /= 0 .and. fj == 0) then
       ff_pair = .true.; fij = fi
    elseif (fi == 0 .and. fj == 0) then
       ff_pair = i < j; fij = 0
    elseif (fi /= 0 .and. fj == -fi) then
       ff_pair = i < j; fij = 0
    endif
  end function ff_pair

  ! initial-state splitting: incoming a (flavour fa_in) emits final-state i
  ! (flavour fi); fa = the flavour of the Born's incoming parton
  logical function if_pair(fa_in, fi, fa)
    integer, intent(in) :: fa_in, fi
    integer, intent(out) :: fa
    if_pair = .true.
    if (fa_in /= 0 .and. fi == 0) then
       fa = fa_in                 ! q -> q + g
    elseif (fa_in /= 0 .and. fi == fa_in) then
       fa = 0                     ! q -> g + q
    elseif (fa_in == 0 .and. fi /= 0) then
       fa = -fi                   ! g -> qbar + q (Born incoming antiquark of q's flavour)
    elseif (fa_in == 0 .and. fi == 0) then
       fa = 0                     ! g -> g + g
    else
       if_pair = .false.; fa = 0
    endif
  end function if_pair

  ! Born: slot 1 = incoming pa, the final-state partons of P4 without slot
  ! 'skip', in order; returns the Born slots of the P4 slots
  subroutine born_slots(P4, fl4, skip, pa, P3, fl3, slot)
    real(dp), intent(in) :: P4(4,8), pa(4)
    integer, intent(in) :: fl4(5), skip
    real(dp), intent(out) :: P3(4,7)
    integer, intent(out) :: fl3(4), slot(5)
    integer :: m, n
    P3 = 0
    P3(:,1) = pa; fl3(1) = fl4(1); slot = 0; slot(1) = 1
    n = 1
    do m = 2, 5
       if (m == skip) cycle
       n = n + 1
       P3(:,n) = P4(:,m); fl3(n) = fl4(m); slot(m) = n
    enddo
    P3(:,5:7) = P4(:,6:8)
  end subroutine born_slots

  subroutine dip_ff(P4, fl4, i, j, k, fij, P3, fl3, val)
    real(dp), intent(in) :: P4(4,8)
    integer, intent(in) :: fl4(5), i, j, k, fij
    real(dp), intent(out) :: P3(4,7), val
    integer, intent(out) :: fl3(4)
    real(dp) :: pij, pik, pjk, y, zi, zj, msq, cc(4,4), mv, cv(4,4), v(4), f
    integer :: slot(5), sij, sk
    pij = mdot(P4(:,i), P4(:,j)); pik = mdot(P4(:,i), P4(:,k)); pjk = mdot(P4(:,j), P4(:,k))
    y = pij/(pij + pik + pjk); zi = pik/(pik + pjk); zj = 1 - zi
    call born_slots(P4, fl4, j, P4(:,1), P3, fl3, slot)
    sij = slot(i); sk = slot(k)
    P3(:,sij) = P4(:,i) + P4(:,j) - y/(1 - y)*P4(:,k)
    P3(:,sk) = P4(:,k)/(1 - y)
    fl3(sij) = fij
    call born31_cc(P3, fl3, msq, cc)
    f = -c16/(2*pij)
    if (fl4(i) /= 0 .and. fl4(j) == 0) then
       val = f*(2/(1 - zi*(1 - y)) - (1 + zi))*cc(sk,sij)
    else
       v = zi*P4(:,i) - zj*P4(:,j)
       call born31_sc(P3, fl3, sij, v, mv, cv)
       if (fl4(i) == 0) then
          val = f*2*((1/(1 - zi*(1 - y)) + 1/(1 - zj*(1 - y)) - 2)*cc(sk,sij) + cv(sk,sij)/pij)
       else
          val = f*TR/CA*(cc(sk,sij) - 2/pij*cv(sk,sij))
       endif
    endif
  end subroutine dip_ff

  subroutine dip_fi(P4, fl4, i, j, fij, P3, fl3, val)
    real(dp), intent(in) :: P4(4,8)
    integer, intent(in) :: fl4(5), i, j, fij
    real(dp), intent(out) :: P3(4,7), val
    integer, intent(out) :: fl3(4)
    real(dp) :: pij, pia, pja, x, zi, zj, msq, cc(4,4), mv, cv(4,4), v(4), f
    integer :: slot(5), sij
    pij = mdot(P4(:,i), P4(:,j)); pia = mdot(P4(:,i), P4(:,1)); pja = mdot(P4(:,j), P4(:,1))
    x = (pia + pja - pij)/(pia + pja); zi = pia/(pia + pja); zj = 1 - zi
    call born_slots(P4, fl4, j, x*P4(:,1), P3, fl3, slot)
    sij = slot(i)
    P3(:,sij) = P4(:,i) + P4(:,j) - (1 - x)*P4(:,1)
    fl3(sij) = fij
    call born31_cc(P3, fl3, msq, cc)
    f = -c16/(2*pij)/x
    if (fl4(i) /= 0 .and. fl4(j) == 0) then
       val = f*(2/(1 - zi + (1 - x)) - (1 + zi))*cc(1,sij)
    else
       v = zi*P4(:,i) - zj*P4(:,j)
       call born31_sc(P3, fl3, sij, v, mv, cv)
       if (fl4(i) == 0) then
          val = f*2*((1/(1 - zi + (1 - x)) + 1/(1 - zj + (1 - x)) - 2)*cc(1,sij) + cv(1,sij)/pij)
       else
          val = f*TR/CA*(cc(1,sij) - 2/pij*cv(1,sij))
       endif
    endif
  end subroutine dip_fi

  subroutine dip_if(P4, fl4, i, k, fa, P3, fl3, val)
    real(dp), intent(in) :: P4(4,8)
    integer, intent(in) :: fl4(5), i, k, fa
    real(dp), intent(out) :: P3(4,7), val
    integer, intent(out) :: fl3(4)
    real(dp) :: pia, pka, pik, x, u, msq, cc(4,4), mv, cv(4,4), v(4), f
    integer :: slot(5), sk
    pia = mdot(P4(:,i), P4(:,1)); pka = mdot(P4(:,k), P4(:,1)); pik = mdot(P4(:,i), P4(:,k))
    x = (pka + pia - pik)/(pka + pia); u = pia/(pia + pka)
    call born_slots(P4, fl4, i, x*P4(:,1), P3, fl3, slot)
    sk = slot(k)
    P3(:,sk) = P4(:,k) + P4(:,i) - (1 - x)*P4(:,1)
    fl3(1) = fa
    call born31_cc(P3, fl3, msq, cc)
    f = -c16/(2*pia)/x
    if (fl4(1) /= 0 .and. fl4(i) == 0) then
       ! q -> q + g
       val = f*(2/(1 - x + u) - (1 + x))*cc(sk,1)
    elseif (fl4(1) == 0 .and. fl4(i) /= 0) then
       ! g -> q(bar) + final q(bar): Born incoming (anti)quark
       val = f*TR/CF*(1 - 2*x*(1 - x))*cc(sk,1)
    else
       v = P4(:,i)/u - P4(:,k)/(1 - u)
       call born31_sc(P3, fl3, 1, v, mv, cv)
       if (fl4(1) /= 0) then
          ! q -> g + final q
          val = f*CF/CA*(x*cc(sk,1) + (1 - x)/x*2*u*(1 - u)/pik*cv(sk,1))
       else
          ! g -> g + g
          val = f*2*((1/(1 - x + u) - 1 + x*(1 - x))*cc(sk,1) + (1 - x)/x*u*(1 - u)/pik*cv(sk,1))
       endif
    endif
  end subroutine dip_if

  pure real(dp) function mdot(a, b)
    real(dp), intent(in) :: a(4), b(4)
    mdot = a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function mdot
end module dip41
