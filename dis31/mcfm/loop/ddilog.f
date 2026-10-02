      function DDILOG(X)
      implicit none
      include 'types.f'
      real(dp):: DDILOG
      integer:: i
      real(dp):: X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      real(dp):: C(0:18),H,ALFA,B0,B1,B2

      DATA ZERO /0.0_dp/, ONE /1.0_dp/
      DATA HALF /0.5_dp/, MALF /-0.5_dp/, MONE /-1.0_dp/, MTWO /-2.0_dp/
      DATA PI3 /3.28986 81336 96453_dp/, PI6 /1.64493 40668 48226_dp/

      DATA C( 0) / 0.42996 69356 08137 0_dp/
      DATA C( 1) / 0.40975 98753 30771 1_dp/
      DATA C( 2) /-0.01858 84366 50146 0_dp/
      DATA C( 3) / 0.00145 75108 40622 7_dp/
      DATA C( 4) /-0.00014 30418 44423 4_dp/
      DATA C( 5) / 0.00001 58841 55418 8_dp/
      DATA C( 6) /-0.00000 19078 49593 9_dp/
      DATA C( 7) / 0.00000 02419 51808 5_dp/
      DATA C( 8) /-0.00000 00319 33412 7_dp/
      DATA C( 9) / 0.00000 00043 45450 6_dp/
      DATA C(10) /-0.00000 00006 05784 8_dp/
      DATA C(11) / 0.00000 00000 86121 0_dp/
      DATA C(12) /-0.00000 00000 12443 3_dp/
      DATA C(13) / 0.00000 00000 01822 6_dp/
      DATA C(14) /-0.00000 00000 00270 1_dp/
      DATA C(15) / 0.00000 00000 00040 4_dp/
      DATA C(16) /-0.00000 00000 00006 1_dp/
      DATA C(17) / 0.00000 00000 00000 9_dp/
      DATA C(18) /-0.00000 00000 00000 1_dp/

      IF(X .EQ. ONE) THEN
       DDILOG=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOG=MALF*PI6
       RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=LOG(-T)
       A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=LOG(-T)
       A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*LOG(T)**2
      END IF

      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END
