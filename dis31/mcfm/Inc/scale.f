
      real(dp):: scale,musq
      common/mcfmscale/scale,musq
!$omp threadprivate(/mcfmscale/)
