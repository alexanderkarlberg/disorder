*
* $Id: dclaus64.F,v 1.2 1996/04/02 16:23:45 mclareni Exp $
*
* $Log: dclaus64.F,v $
* Revision 1.2  1996/04/02 16:23:45  mclareni
* More precise dclaus64 (C326), test added and C344 removed from TESTALL
*
* Revision 1.1.1.1  1996/04/01 15:02:03  mclareni
* Mathlib gen
*
*
      function DCLAUS(X)
      implicit none
      include 'types.f'
      real(dp):: DCLAUS,H,V,U,ALFA,B0,B1,B2,S,X
      integer:: i
      real(dp):: A(0:8),B(0:13),R1,HF,PI,PI2,PIH,RPIH
 
      PARAMETER (R1 = 1._dp, HF =R1/2._dp)
      PARAMETER (PI = 3.14159 26535 89793 24_dp)
      PARAMETER (PI2 = 2._dp*PI, PIH = PI/2._dp, RPIH = 2._dp/PI)
 
      DATA A( 0) / 0.02795 28319 73575 6613_dp/
      DATA A( 1) / 0.00017 63088 74389 8116_dp/
      DATA A( 2) / 0.00000 12662 74146 1157_dp/
      DATA A( 3) / 0.00000 00117 17181 8134_dp/
      DATA A( 4) / 0.00000 00001 23006 4129_dp/
      DATA A( 5) / 0.00000 00000 01395 2729_dp/
      DATA A( 6) / 0.00000 00000 00016 6908_dp/
      DATA A( 7) / 0.00000 00000 00000 2076_dp/
      DATA A( 8) / 0.00000 00000 00000 0027_dp/
 
      DATA B( 0) / 0.63909 70888 57265 341_dp/
      DATA B( 1) /-0.05498 05693 01851 716_dp/
      DATA B( 2) /-0.00096 12619 45950 606_dp/
      DATA B( 3) /-0.00003 20546 86822 550_dp/
      DATA B( 4) /-0.00000 13294 61695 426_dp/
      DATA B( 5) /-0.00000 00620 93601 824_dp/
      DATA B( 6) /-0.00000 00031 29600 656_dp/
      DATA B( 7) /-0.00000 00001 66351 954_dp/
      DATA B( 8) /-0.00000 00000 09196 527_dp/
      DATA B( 9) /-0.00000 00000 00524 004_dp/
      DATA B(10) /-0.00000 00000 00030 580_dp/
      DATA B(11) /-0.00000 00000 00001 820_dp/
      DATA B(12) /-0.00000 00000 00000 110_dp/
      DATA B(13) /-0.00000 00000 00000 007_dp/
 
      V=MOD(ABS(X),PI2)
      S=SIGN(R1,X)
      IF(V .GT. PI) THEN
       V=PI2-V
       S=-S
      ENDIF
      IF(V .EQ. 0._dp .OR. V .EQ. PI) THEN
       H=0._dp
      ELSEIF(V .LT. PIH) THEN
       U=RPIH*V
       H=2._dp*U**2-1._dp
       ALFA=H+H
       B1=0._dp
       B2=0._dp
       DO 1 I = 8,0,-1
       B0=A(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=V*(1._dp-LOG(V)+HF*V**2*(B0-H*B2))
      ELSE
       U=RPIH*V-2._dp
       H=2._dp*U**2-1._dp
       ALFA=H+H
       B1=0._dp
       B2=0._dp
       DO 2 I = 13,0,-1
       B0=B(I)+ALFA*B1-B2
       B2=B1
    2  B1=B0
       H=(PI-V)*(B0-H*B2)
      ENDIF
      DCLAUS=S*H
      RETURN
      END
