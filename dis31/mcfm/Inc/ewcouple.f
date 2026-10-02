
      real(dp):: Gf,gw,xw,gwsq,esq,vevsq
      common/ewcouple/Gf,gw,xw,gwsq,esq,vevsq
!$omp threadprivate(/ewcouple/)
