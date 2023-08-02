! Program to merge grids after stage 1
      program mergegrids
      implicit none
      integer maxfiles
      parameter (maxfiles=9999)
      character *(100) files(maxfiles)
      integer ifile,nfiles,i,j
      integer ndim, ndmx, mxdim
      parameter (ndim = 7, ndmx=50, mxdim=7)
      double precision jj, xin(ndmx,mxdim,maxfiles), xi(ndmx,mxdim),
     $     xval
      
      nfiles = -1
      
      do ifile=1,maxfiles
         CALL getarg(ifile,files(ifile))         
         if(trim(files(ifile)).eq.'') then 
            nfiles=ifile-1
            write(6,*) 'mergedata.exe found',nfiles,
     $           'files on the command line ...'
            exit
         endif
      enddo
      if(nfiles.eq.-1) stop "Increase maxfiles"

      do ifile=1,nfiles
         open(unit=11,file=files(ifile),status='old')
         write(6,*)'* Reading in grid from ',trim(files(ifile)),' *'
         call flush(6)
         do j=1,ndim
            read(11,203) jj,(xin(i,j,ifile),i=1,ndmx)
         enddo
         close(11)
      enddo
      xi = 0d0
      do ifile=1,nfiles
         do j=1,ndim
            do i=1,ndmx
               xi(i,j) = xi(i,j) + xin(i,j,ifile)
            enddo
         enddo
      enddo
      xi = xi/nfiles

      open(unit=11,file='merged-grids.dat',status='unknown')
      open(unit=111,file='merged-grids.top',status='unknown')
      write(6,*)'* Writing out grid to ','merged-grids.dat',' *'
      call flush(6)
      do j=1,ndim
         write(11,203) jj,(xi(i,j),i=1,ndmx)
         write(111,*) '# XDIM index  ', j-1
         do i=1,ndmx
            xval=(i*1d0)/ndmx
            write(111,'(D16.8,D16.8)') xval, xi(i,j)
         enddo
         write(111,*) ''
         write(111,*) ''
      enddo
      close(11)
      close(111)

      
 203  FORMAT(/(5z16))
      end





