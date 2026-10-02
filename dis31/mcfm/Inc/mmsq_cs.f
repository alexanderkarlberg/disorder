
      real(dp):: mmsq_cs(0:2,2,2)
      common/mmsq_cs/mmsq_cs
!$omp threadprivate(/mmsq_cs/)
