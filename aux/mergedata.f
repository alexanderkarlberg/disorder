c program to merge gnuplot data tables.
c a line starting with '#' followed by a line with 4 numbers
c is considered the beginning of a data set. All files to be
c merged must be identical in format.
c When the program starts, it expects as input a list of filenames,
c one per line, to be merged. An empty lines terminates the list.
      implicit none
      integer maxfiles,maxlines
      parameter (maxfiles=1000,maxlines=25000)
      character *(100) files(maxfiles)
      character *(500) line(maxlines,maxfiles)
      integer nlines(maxfiles)
      integer ifile,nfiles,ios,k,imethod
      character cmethod
      real * 8 v1,v2,v3,v4,y,err
      integer ilength
      external ilength

C - kh modification started here >>>>>>>
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
      elseif(cmethod.eq.' ') then
         imethod=0 ! If nothing follows ./mergedata.exe on the command
                   ! line this setting is acquired which then steers the
                   ! code to read in the combination mode and input files
                   ! from the read command prompts as it had been doing.
      else
         write(6,*) 'Combination mode must be "1-5" or " " : ',
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
     $                    'files on the command line ...'
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
      read(*,*) imethod
      write(*,*) ' enter files'
      do ifile=1,maxfiles
         read(*,'(a)') files(ifile)
         if(files(ifile).eq.' ') then
            nfiles=ifile-1
            goto 10
         endif
      enddo
      write(*,*) ' too manny files, increase maxfiles'
      call exit(-1)
      endif
 10   continue
c load data
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
      do k=1,nlines(1)
         read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3,v4
         if(ios.ne.0) then
            write(12,'(a)') line(k,1)(1:ilength(line(k,1)))
         else
            if(imethod.eq.1) then
               y=v3
               err=v4**2
            elseif(imethod.eq.2) then
               if(v4.ne.0) then
                  y=v3/v4**2
                  err=1/v4**2
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
            endif
            do ifile=2,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               if(imethod.eq.1.or.imethod.eq.3) then
                  y=y+v3
                  err=err+v4**2
               elseif(imethod.eq.2) then
                  if(v4.ne.0) then
                     y=y+v3/v4**2
                     err=err+1/v4**2
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
               endif
            enddo
            if(imethod.eq.1) then
               y=y/nfiles
               err=sqrt(err/nfiles**2)
            elseif(imethod.eq.2) then
               if(err.ne.0) then
                  y=y/err
                  err=1/sqrt(err)
               else
                  y=0
                  err=0
               endif
            elseif(imethod.ge.3) then
               err=sqrt(err)
            endif
            write(12,'(4(1x,d14.8))') v1,v2,y,err
         endif
      enddo
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
