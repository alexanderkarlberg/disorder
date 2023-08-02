      program mergedata

c program to merge gnuplot data tables.
c a line starting with '#' followed by a line with 4 numbers
c is considered the beginning of a data set. All files to be
c merged must be identical in format.
c When the program starts, it expects as input a list of filenames,
c one per line, to be merged. An empty lines terminates the list.
      implicit none
      integer maxfiles,maxlines
      parameter (maxfiles=1000,maxlines=10000)
      character *(100) files(maxfiles)
      character *(100) line(maxlines,maxfiles)
      integer nlines(maxfiles), int_fraction
      integer ifile,nfiles,ios,k,imethod,N_NaN,npoints,ios2
      character *(2) cmethod
      real * 8 v1,v2,v3,v4,y,err,integral,xsec,scale, fraction
      real * 8 central, error, outlier, db_npoints, db_nfiles
      real * 8  sqrt_nfiles, xmin, xmax, xbin, var
      real * 8 data_val(maxfiles,maxlines),error_val(maxfiles,maxlines)
      real * 8 median, IQR, limit
      integer sorted_index(maxfiles,maxlines), first_index, last_index
      integer ilength,indexx,i,j
      external ilength, integral,integral_error
      character *(100), c1, c2, c3
      integer i4, iindex, ibin, nbins
      parameter (nbins=100)
      real * 8  hbin(nbins)
      real * 8 vp1,vp2,vp3,vp4,integral_error
      
      fraction=0d0
      int_fraction=0
      outlier = 100d0 ! Number of STD away from central value to be considered an outlier 
C     - kh modification started here >>>>>>>
      imethod=-1
      CALL getarg(1,cmethod)
      if(cmethod.eq.'1') then
         imethod=1
      elseif(cmethod.eq.'2') then
         imethod=2
      elseif(cmethod.eq.'3') then
         imethod=3
      elseif(cmethod.eq.'4') then
         imethod=4
      elseif(cmethod.eq.'5') then
         imethod=5
      elseif(cmethod.eq.'6') then
         imethod=6
      elseif(cmethod.eq.'7') then
         imethod=7
      elseif(cmethod.eq.'8') then
         imethod=8
      elseif(cmethod.eq.'9') then
         imethod=9
      elseif(cmethod.eq.'10') then
         imethod=10
      elseif(cmethod.eq.'11') then
         imethod=11
      elseif(cmethod.eq.'12') then
         imethod=12
      elseif(cmethod.eq.'13') then
         imethod=13
      elseif(cmethod.eq.'14') then
         imethod=14
      elseif(cmethod.eq.'15') then
         imethod=15
      elseif(cmethod.eq.'16') then
         imethod=16
      elseif(cmethod.eq.'17') then
         imethod=17
      elseif(cmethod.eq.' ') then
         imethod=0 ! If nothing follows ./mergedata.exe on the command
                   ! line this setting is acquired which then steers the
                   ! code to read in the combination mode and input files
                   ! from the read command prompts as it had been doing.
      else
         write(6,*) 'Combination mode must be "1-17" or " " : ',
     $               cmethod
         write(6,*) 'Quitting ...'
         stop
      endif
      if(imethod.ne.0) then
         do ifile=1,maxfiles
            CALL getarg(ifile+1,files(ifile))         
            if(trim(files(ifile)).eq.'') then 
               nfiles=ifile-1
               write(6,*) 'mergedata.exe found',nfiles,
     $              'files on the command line ...'
               goto 9
            endif
         enddo
      endif
 9    continue

      if(imethod.eq.0) then
C - <<<<<<< kh modification ended here.
      write(*,*) ' enter 1 for combining sets with equal statistics'
      write(*,*) ' 2 to combine uneven sets'
      write(*,*) ' 3 to add sets (like born+virtual+real ... etc'
      write(*,*) ' 4 to get maximum'
      write(*,*) ' 5 to get minimum'
      write(*,*) ' 6 to normalise all distributions'
      write(*,*) ' 7 to accumulate all distributions'
      write(*,*) ' 8 to multiply with scalar (asks for input)'
      write(*,*) ' 9 to divide with scalar (asks for input)'
      write(*,*) ' 10 to divide two sets'
      write(*,*) ' 11 to multiply two sets'
      write(*,*) ' 12 to combine uneven sets discarding outliers'
      write(*,*) ' 13 compute the truncated mean of the 
     $     (2%-98%) middle sample'
      write(*,*) ' 14 to do weighted average with abs(error)'
      write(*,*) ' 15 to do median after removal of outliers'
      write(*,*) ' 16 to do mean after removal of outliers'
      write(*,*) ' 17 to accumulate all distributions starting 
     $     from the end'
      read(*,*) imethod
      write(*,*) ' enter files'
      do ifile=1,maxfiles
         read(*,'(a)') files(ifile)
         if(files(ifile).eq.' ') then
            nfiles=ifile-1
            if(imethod.eq.6.and.nfiles.ne.1) then
               print*, 'Called with',nfiles, 'files and method', imethod
               print*, 'Should be called with only 1 file'
               print*, 'Exiting...'
               stop
            endif
            if((imethod.eq.8).or.(imethod.eq.9).and.nfiles.ne.1) then
               print*, 'Called with',nfiles, 'files and method', imethod
               print*, 'Should be called with 1 file and 1 scalar'
               print*, 'Exiting...'
               stop
            endif
            goto 10
         endif
      enddo
      write(*,*) ' too many files, increase maxfiles'
      call exit(-1)
      endif
 10   continue
c load data
      db_nfiles=nfiles
      sqrt_nfiles=sqrt(db_nfiles)
      do ifile=1,nfiles
         open(unit=11,file=files(ifile),status='old')
         do k=1,maxlines+1
            read(unit=11,fmt='(a)',end=111) line(k,ifile)
            if(k.eq.maxlines+1) then
               write(*,*) ' too many lines in file, increase maxlines'
               call exit(-1)
            endif
            goto 12
 111        nlines(ifile)=k-1
            goto 11
 12         continue
         enddo
 11      continue
      enddo
      do ifile=1,nfiles
         if(nlines(ifile).ne.nlines(1)) then
            write(*,*) ' error: file', files(ifile),
     1           ' does not match in length'
            call exit(-1)
         endif
      enddo
      if((imethod.eq.8).or.(imethod.eq.9)) then
         write(*,*) 'Input scale to multiply/divide by:'
         read(*,*) scale
         write(*,*), 'Scale is:', scale
      endif
      if(imethod.eq.12) then
         write(*,*) 'Input tolerance for outliers 
     $        (small number throws away more outliers):'
         read(*,*) outlier
         write(*,*) 'Tolerance is:', outlier
      endif
      if(imethod.eq.13) then
c     Fill arrays with data for sorting
         do k=1,nlines(1)
            do ifile=1,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               if(ios.eq.0) then
                  data_val(ifile,k) = v3
                  error_val(ifile,k) = v4
               endif
            enddo
         enddo

c     Sort data
         do k=1,nlines(1)
            call sortrx(nfiles,data_val(:,k),sorted_index(:,k))
         enddo
         do i=1,100
            limit=(0.5d0*i)/100d0
         do k=1,nlines(1)
c            call sortrx(nfiles,data_val(:,k),sorted_index(:,k))
c     Compute median
            if(mod(nfiles,2).eq.0) then !nfiles even
               median = (data_val(sorted_index(nfiles/2,k),k)+
     $              data_val(sorted_index(nfiles/2+1,k),k))/2d0
            else
               median = data_val(sorted_index(nfiles/2+1,k),k)
            endif
            IQR = data_val(sorted_index((3*nfiles)/4,k),k)-
     $           data_val(sorted_index((nfiles)/4,k),k)
            read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
            if(ios.ne.0) then
               write(i+100,'(a)') line(k,1)(1:ilength(line(k,1)))
c               write(120,'(a)') line(k,1)(1:ilength(line(k,1)))
            else
               y=0d0
               err=0d0
               npoints=0
               read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
               do ifile=1,nfiles
                  if((1d0*ifile)/(1d0*nfiles)
     $                 .gt.(0.5d0-1d0*limit)
     $                 .and.
     $                 (1d0*ifile)/(1d0*nfiles)
     $                 .le.(0.5d0+1d0*limit)) then ! Take only middle $limit %
                     npoints=npoints+1
                     y=y+data_val(sorted_index(ifile,k),k)
                     err=err+error_val(sorted_index(ifile,k),k)**2
                  endif
               enddo
               y=y/npoints
               err=sqrt(err)/npoints
               write(i+100,'(5(1x,d16.8))') v1,v2,y,err,1d0*limit
c               write(120,'(4(1x,d16.8))') v1,v2,median,IQR
            endif
         enddo
         enddo
         stop
      endif
      if(imethod.eq.15) then
c     Fill arrays with data for sorting
c         print*, 'Filling arrays'
         do k=1,nlines(1)
            do ifile=1,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               if(ios.eq.0) then
                  data_val(ifile,k) = v3
                  error_val(ifile,k) = v4
               endif
            enddo
         enddo

         do k=1,nlines(1)
         first_index = 0
         last_index = 0
c     Sort data
            call sortrx(nfiles,data_val(:,k),sorted_index(:,k))
c     Compute median
            if(mod(nfiles,2).eq.0) then !nfiles even
               median = (data_val(sorted_index(nfiles/2,k),k)+
     $              data_val(sorted_index(nfiles/2+1,k),k))/2d0
            else
               median = data_val(sorted_index(nfiles/2+1,k),k)
            endif
c     Compute Inter Quartile Range (central 50%)
            IQR = data_val(sorted_index((3*nfiles)/4,k),k)-
     $           data_val(sorted_index((nfiles)/4,k),k)
c     Find first index outside 3 times half IQR and last index inside.
            do i = 1, nfiles
c               print*, 'In loop doing i =', i
               if(data_val(sorted_index(i,k),k).lt.
     $              (median - (IQR/2d0)*10d0)) then
                  first_index = first_index + 1
               endif
               if(data_val(sorted_index(i,k),k).lt.
     $              (median + (IQR/2d0)*10d0)) then
                  last_index = last_index + 1
               endif
            enddo
c            print*, 'median, IQR/2d0', median, IQR/2d0
c            print*, 'first_index', first_index
c            print*, 'last_index', last_index
c            print*, 'nfiles', nfiles
c            stop
c     Recompute the median having thrown out outliers
            if(mod(last_index-first_index,2).eq.0) then !nfiles even
               median = (data_val(sorted_index(
     $              (last_index-first_index)/2,k),k)+
     $              data_val(sorted_index(
     $              (last_index-first_index)/2+1,k),k))/2d0
            else
               median = data_val(sorted_index(
     $              (last_index-first_index)/2+1,k),k)
            endif
            
            read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
            if(ios.ne.0) then
               write(12,'(a)') line(k,1)(1:ilength(line(k,1)))
            else
               y=0d0
               err=0d0
               npoints=0
               read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
               write(12,'(4(1x,d16.8))') v1,v2,median,IQR
     $              /sqrt(nfiles*1d0)
            endif
         enddo
         return
      endif

      if(imethod.eq.16) then
c     Fill arrays with data for sorting
c         print*, 'Filling arrays'
         do k=1,nlines(1)
            do ifile=1,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               if(ios.eq.0) then
                  data_val(ifile,k) = v3
                  error_val(ifile,k) = v4
               endif
            enddo
         enddo

c     Sort data
         do k=1,nlines(1)
            call sortrx(nfiles,data_val(:,k),sorted_index(:,k))
         enddo
         do i =1,100
            limit=(10d0*i)/10d0
            do k=1,nlines(1)
               first_index = 0
               last_index = 0
c     Compute median
               if(mod(nfiles,2).eq.0) then !nfiles even
                  median = (data_val(sorted_index(nfiles/2,k),k)+
     $                 data_val(sorted_index(nfiles/2+1,k),k))/2d0
               else
                  median = data_val(sorted_index(nfiles/2+1,k),k)
               endif
c     Compute Inter Quartile Range (central 68%)
               IQR = data_val(sorted_index((84*nfiles)/100,k),k)-
     $              data_val(sorted_index((16*nfiles)/100,k),k)
c     Find first index outside 3 times half IQR and last index inside.
               do ifile = 1, nfiles
c     print*, 'In loop doing i =', i
                  if(data_val(sorted_index(ifile,k),k).lt.
     $                 (median - (IQR/2d0)*limit)) then
                     first_index = first_index + 1
                  endif
                  if(data_val(sorted_index(ifile,k),k).lt.
     $                 (median + (IQR/2d0)*limit)) then
                     last_index = last_index + 1
                  endif
               enddo
c     Recompute the median having thrown out outliers
               
               read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
               if(ios.ne.0) then
                  write(i+100,'(a)') line(k,1)(1:ilength(line(k,1)))
               else
                  y = 0d0
                  err = 0d0
                  do j = first_index, last_index
                     y = y + data_val(sorted_index(j,k),k)
                     err = err + error_val(sorted_index(j,k),k)**2
                  enddo
                  y = y/(last_index-first_index)
                  err=sqrt(err)/(last_index-first_index) 
c     y=sum(data_val(first_index:last_index,k),dim=1)/
c     $              (last_index-first_index)
                  if(y.ne.y) y = 0d0
c     err=sum(data_val(first_index:last_index,k)**2,dim=1)
c     err=sqrt(err/(last_index-first_index))
                  if(err.ne.err) err = 0d0
                  npoints=0
                  read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
                  write(i+100,'(6(1x,d16.8))') v1,v2,y,err, limit*1d0, 
     $                 (1d0*last_index-1d0*first_index)/(1d0*nfiles)
               endif
            enddo
         enddo
            return
         endif

      do k=1,nlines(1)
         N_NaN=0
         read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
         read(unit=line(k-1,1),fmt=*,iostat=ios2) vp1,vp2,vp3,vp4
         if(ios.ne.0) then
            y=0d0
            err=0d0
            write(12,'(a)') line(k,1)(1:ilength(line(k,1)))
            if(imethod.eq.6.or.imethod.eq.17) then
               read(unit=line(k,1),fmt=*,iostat=ios) c1,c2,c3,i4
               xsec=integral(line(k:,1),nlines(1)) !Compute integral of distribution
               err=integral_error(line(k:,1),nlines(1)) !Compute integral of error
               y=xsec
            endif
         else
            if(imethod.eq.1) then
               if(v3.ne.v3.or.v4.ne.v4.or.abs(v3).ge.1d100) then !Tests for NaN and infinities
                  write(*,*) 'Found NaN in line', k, 'of first file'
                  N_NaN=N_NaN+1
               else
                  y=v3
                  err=v4**2
               endif
            elseif(imethod.eq.2) then
               if(v4.ne.0) then
                  y=v3/v4**2
                  err=1/v4**2
               else
                  y=0
                  err=0
               endif
            elseif(imethod.eq.14) then
               if(v4.ne.0) then
                  y=v3/abs(v4)
                  err=1/abs(v4)
               else
                  y=0
                  err=0
               endif
            elseif(imethod.eq.3) then
               y=v3
               err=v4**2
            elseif(imethod.eq.4) then
               y=v3
               err=v4**2               
            elseif(imethod.eq.5) then
               y=v3
               err=v4**2
            elseif(imethod.eq.6) then
               y=v3/xsec
               err=(v4)**2/xsec
            elseif(imethod.eq.7) then
               y=v3+y
               err=err**2+(v4)**2
            elseif(imethod.eq.17) then
               if(ios2.eq.0) then
                  y=y-vp3*(vp2-vp1)
                  err=err**2-(vp4*(vp2-vp1))**2
!                  err=(err-vp4*(vp2-vp1))**2
               else
                  y = y
                  err = err**2
               endif
            elseif(imethod.eq.8) then
               y=v3*scale
               err=(scale*v4)**2
            elseif(imethod.eq.9) then
               y=v3/scale
               err=(v4/scale)**2
            elseif(imethod.eq.10) then
                  y=v3
                  err=(v4/v3)**2
            elseif(imethod.eq.11) then
                  y=v3
                  err=(v4/v3)**2
            elseif(imethod.eq.12) then
               if(v3.ne.v3.or.v4.ne.v4.or.abs(v3).ge.1d100) then !Tests for NaN and infinities
                  write(*,*) 'Found NaN in line', k, 'of first file'
                  N_NaN=N_NaN+1
               else
                  y=v3
                  err=v4**2
               endif
            endif
            do ifile=2,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               if(imethod.eq.1.or.imethod.eq.3) then
                  if(v3.ne.v3.or.v4.ne.v4.or.abs(v3).ge.1d100) then !Tests for NaN and infinities
                     write(*,*) 'Found NaN in line', k, 'of file', ifile
                     N_NaN=N_NaN+1
                  else
                     y=y+v3
                     err=err+v4**2
                  endif
               elseif(imethod.eq.2) then
                  if(v4.ne.0) then
                     y=y+v3/v4**2
                     err=err+1/v4**2
                  endif
               elseif(imethod.eq.14) then
                  if(v4.ne.0) then
                     y=y+v3/abs(v4)
                     err=err+1/abs(v4)
                  endif
               elseif(imethod.eq.3) then
                  y=y+v3
                  err=err+v4**2
               elseif(imethod.eq.4) then
                  if(v3.gt.y) then
                     y=v3
                     err=v4**2
                  endif
               elseif(imethod.eq.5) then
                  if(v3.lt.y) then
                     y=v3
                     err=v4**2
                  endif
               elseif(imethod.eq.10) then
                  if(v3.ne.0d0) then
                     y=y/v3
                     err=y**2*(err+(v4/v3)**2)
                  else
                     y=0d0
                     err=0d0
                  endif
               elseif(imethod.eq.11) then
                  y=y*v3
                  err=y**2*(err+(v4/v3)**2)
               elseif(imethod.eq.12) then
                  if(v3.ne.v3.or.v4.ne.v4.or.abs(v3).ge.1d100) then !Tests for NaN and infinities
                     write(*,*) 'Found NaN in line', k, 'of file', ifile
                     N_NaN=N_NaN+1
                  else
                     y=y+v3
                     err=err+v4**2
                  endif
               endif
            enddo
            if(imethod.eq.1.and.(nfiles.ne.N_NaN).or.
     $           imethod.eq.12.and.(nfiles.ne.N_NaN)) then
               y=y/(nfiles-N_NaN)
               err=sqrt(err/(nfiles-N_Nan)**2)
            elseif(imethod.eq.2) then
               if(err.ne.0) then
                  y=y/err
                  err=1/sqrt(err)
               else
                  y=0d0
                  err=0d0
               endif
            elseif(imethod.eq.14) then
               if(err.ne.0) then
                  y=y/err
                  err=1/err
               else
                  y=0d0
                  err=0d0
               endif
c     elseif(imethod.ge.3) then
            else
               err=sqrt(abs(err))
            endif
            
            if(imethod.eq.12) then
               central=y
               error=err
               y=0d0
               err=0d0
               npoints=0
               do ifile=1,nfiles
                  read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
                  if(v3.ne.v3.or.v4.ne.v4.or.abs(v3).ge.1d100) then !Tests for NaN and infinities
                     write(*,*) 'Found NaN in line', k, 'of file', ifile
                     N_NaN=N_NaN+1
                  elseif((abs((v3-central)/(error
     $                    *sqrt_nfiles)).lt.outlier).and.
     $                    v4.lt.outlier*error*sqrt_nfiles) then
                     y=y+v3
                     err=err+v4**2
                     npoints=npoints+1
                  endif
               enddo
               db_npoints=npoints
               if(y.ne.0d0) then
c                  print*, 'Ratio of points discarded', 
c     $                 (db_nfiles-db_npoints)/db_nfiles
                  fraction = (db_nfiles-db_npoints)/db_nfiles + fraction
                  int_fraction=int_fraction+1
               endif
               if(npoints.gt.0) then
                  y=y/(npoints)
                  err=sqrt(err/(npoints)**2)
               endif
            endif
            write(12,'(4(1x,d16.8))') v1,v2,y,err
         endif
      enddo
      if(imethod.eq.12) then
         print*, 'Average number of outliers:', 
     $        fraction/int_fraction*100d0, '%'
      endif
      end


      function ilength(line)
      integer ilength
      character *(*) line
      ilength=len(line)
      do j=ilength,1,-1
         if(line(j:j).ne.' ') then
            ilength=j
            return
         endif
      enddo
      ilength=0
      end

      function integral(line, nlines)
      implicit none
      real * 8 integral
      integer maxfiles,maxlines
      parameter (maxfiles=10000,maxlines=1500)
      character *(100) file
      character *(100) line(maxlines)
      integer nlines,index, indexx
      integer ios,k
      real * 8 v1,v2,v3,v4
      character *(100), c1, c2, c3
      integer i4
      integer ilength
      external ilength
      index=0
      
      integral = 0d0
      do k=2,nlines
         read(unit=line(k),fmt=*,iostat=ios) v1,v2,v3,v4
         if(ios.ne.0) then
            return
         else
            integral=integral+v3*(v2-v1)
         endif
      enddo

      end

      function integral_error(line, nlines)
      implicit none
      real * 8 integral_error
      integer maxfiles,maxlines
      parameter (maxfiles=10000,maxlines=1500)
      character *(100) file
      character *(100) line(maxlines)
      integer nlines,index, indexx
      integer ios,k
      real * 8 v1,v2,v3,v4
      character *(100), c1, c2, c3
      integer i4
      integer ilength
      external ilength
      index=0
      
      integral_error = 0d0
      do k=2,nlines
         read(unit=line(k),fmt=*,iostat=ios) v1,v2,v3,v4
         if(ios.ne.0) then
            integral_error = sqrt(integral_error)
            return
         else
            integral_error=integral_error+(v4*(v2-v1))**2
         endif
      enddo


      end


C From Leonard J. Moss of SLAC:

C Here's a hybrid QuickSort I wrote a number of years ago.  It's
C based on suggestions in Knuth, Volume 3, and performs much better
C than a pure QuickSort on short or partially ordered input arrays.  

      SUBROUTINE SORTRX(N,DATA,INDEX)
C===================================================================
C
C     SORTRX -- SORT, Real input, indeX output
C
C
C     Input:  N     INTEGER
C             DATA  REAL
C
C     Output: INDEX INTEGER (DIMENSION N)
C
C This routine performs an in-memory sort of the first N elements of
C array DATA, returning into array INDEX the indices of elements of
C DATA arranged in ascending order.  Thus,
C
C    DATA(INDEX(1)) will be the smallest number in array DATA;
C    DATA(INDEX(N)) will be the largest number in DATA.
C
C The original data is not physically rearranged.  The original order
C of equal input values is not necessarily preserved.
C
C===================================================================
C
C SORTRX uses a hybrid QuickSort algorithm, based on several
C suggestions in Knuth, Volume 3, Section 5.2.2.  In particular, the
C "pivot key" [my term] for dividing each subsequence is chosen to be
C the median of the first, last, and middle values of the subsequence;
C and the QuickSort is cut off when a subsequence has 9 or fewer
C elements, and a straight insertion sort of the entire array is done
C at the end.  The result is comparable to a pure insertion sort for
C very short arrays, and very fast for very large arrays (of order 12
C micro-sec/element on the 3081K for arrays of 10K elements).  It is
C also not subject to the poor performance of the pure QuickSort on
C partially ordered data.
C
C Created:  15 Jul 1986  Len Moss
C
C===================================================================
 
      INTEGER   N,INDEX(N)
      REAL * 8     DATA(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL * 8     DATAP
 
C     QuickSort Cutoff
C
C     Quit QuickSort-ing when a subsequence contains M or fewer
C     elements and finish off at end with straight insertion sort.
C     According to Knuth, V.3, the optimum value of M is around 9.
 
      INTEGER   M
      PARAMETER (M=9)
 
C===================================================================
C
C     Make initial guess for INDEX
 
      DO 50 I=1,N
         INDEX(I)=I
   50    CONTINUE
 
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
 
      IF (N.LE.M) GOTO 900
 
C===================================================================
C
C     QuickSort
C
C     The "Qn:"s correspond roughly to steps in Algorithm Q,
C     Knuth, V.3, PP.116-117, modified to select the median
C     of the first, last, and middle elements as the "pivot
C     key" (in Knuth's notation, "K").  Also modified to leave
C     data in place and produce an INDEX array.  To simplify
C     comments, let DATA[I]=DATA(INDEX(I)).
 
C Q1: Initialize
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
C Q2: Sort the subsequence DATA[L]..DATA[R].
C
C     At this point, DATA[l] <= DATA[m] <= DATA[r] for all l < L,
C     r > R, and L <= m <= R.  (First time through, there is no
C     DATA for l < L or r > R.)
 
      I=L
      J=R
 
C Q2.5: Select pivot key
C
C     Let the pivot, P, be the midpoint of this subsequence,
C     P=(L+R)/2; then rearrange INDEX(L), INDEX(P), and INDEX(R)
C     so the corresponding DATA values are in increasing order.
C     The pivot key, DATAP, is then DATA[P].
 
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
 
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
 
C     Now we swap values between the right and left sides and/or
C     move DATAP until all smaller values are on the left and all
C     larger values are on the right.  Neither the left or right
C     side will be internally ordered yet; however, DATAP will be
C     in its final position.
 
  300 CONTINUE
 
C Q3: Search for datum on left >= DATAP
C
C     At this point, DATA[L] <= DATAP.  We can therefore start scanning
C     up from L, looking for a value >= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         I=I+1
         IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
C Q4: Search for datum on right <= DATAP
C
C     At this point, DATA[R] >= DATAP.  We can therefore start scanning
C     down from R, looking for a value <= DATAP (this scan is guaranteed
C     to terminate since we initially placed DATAP near the middle of
C     the subsequence).
 
         J=J-1
         IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
 
C Q5: Have the two scans collided?
 
      IF (I.LT.J) THEN
 
C Q6: No, interchange DATA[I] <--> DATA[J] and continue
 
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
 
C Q7: Yes, select next subsequence to sort
C
C     At this point, I >= J and DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     for all L <= l < I and J < r <= R.  If both subsequences are
C     more than M elements long, push the longer one on the stack and
C     go back to QuickSort the shorter; if only one is more than M
C     elements long, go back and QuickSort it; otherwise, pop a
C     subsequence off the stack and QuickSort it.
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
C===================================================================
C
C Q9: Straight Insertion sort
 
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX(P+1) = INDEX(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
C===================================================================
C
C     All done
 
      END
