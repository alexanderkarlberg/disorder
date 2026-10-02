
      real(dp) :: md,mu,ms,mc,mb,mt
      real(dp) :: mel,mmu,mtau
      real(dp) :: hmass,hwidth
      real(dp) :: wmass,wwidth
      real(dp) :: zmass,zwidth
      real(dp) :: twidth
      real(dp) :: tauwidth
      real(dp) :: mtausq,mcsq,mbsq
      common/masses/md,mu,ms,mc,mb,mt,mel,mmu,mtau,hmass,hwidth,wmass,wwidth,zmass,zwidth,twidth,tauwidth,mtausq,mcsq,mbsq
!$omp threadprivate(/masses/)
