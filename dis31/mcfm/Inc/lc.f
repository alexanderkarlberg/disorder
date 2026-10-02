
      integer:: colourchoice
!--- 'colourchoice' allows calculation by colour structure
!--- For Gflag=.true. [QQGG, QQGGG processes]
!--- 1) Only leading colour ( NCF . N )
!--- 2) Only sub-leading ( NCF . 1/N )
!--- 3) Only sub-sub-leading ( NCF . 1/N**3 )  [QQGGG only]
!--- 0) The total
!--- For Qflag=.true. [QQBQQB process]
!--- 1) Only leading colour ( NCF . 1 )
!--- 2) Only sub-leading ( NCF . 1/N ) [Identical quarks only]
!--- 0) The total
      common/ColC/colourchoice
!$omp threadprivate(/ColC/)
