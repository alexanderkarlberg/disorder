
      real(dp):: s(mxpart,mxpart)
      common/sprods/s
!$omp threadprivate(/sprods/)
