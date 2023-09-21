c program to compute pdfuncertainties from a list of datafiles a line
c starting with '#' followed by a line with 4 numbers is considered the
c beginning of a data set. The file returned has the central PDF value
c in column 3 and the symmetric PDF error in column 4.
      implicit none
      integer maxfiles,maxlines
      parameter (maxfiles=1000,maxlines=2500)
      character *(100) files(maxfiles)
      character *(500) line(maxlines,maxfiles)
      integer nlines(maxfiles),nset
      integer ifile,nfiles,ios,k,ipdf,nmempdf_end
      character *(100) cpdf
      real * 8 v1,v2,v3,v4
      double precision res(0:200), central, errminus, errplus,
     $     errsymm
      integer ilength
      external ilength


      CALL getarg(1,cpdf)

      if(cpdf.eq.''.or.cpdf.eq.'--help') then
C - <<<<<<< kh modification ended here.
         write(*,*) ' Syntax ./getpdfuncert PDFNUMBER list_of_files'
         write(*,*)
     $    ' Example: ./getpdfuncert 91200 data1.dat data2.dat data3.dat'

      endif

      call initPDFSetByName(trim(cpdf))
      call numberPDF(nmempdf_end)

      OPEN(UNIT=51, FILE='pdfuncert.dat', ACTION="write")

      res = 0d0
      
      do ifile=1,maxfiles
         CALL getarg(ifile+1,files(ifile))         
         if(trim(files(ifile)).eq.'') then 
            nfiles=ifile-1
            write(6,*) 'mergedata.exe found',nfiles,
     $           'files on the command line ...'
            goto 9
         endif
      enddo
 9    continue

      
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
            write(51,'(a)') line(k,1)(1:ilength(line(k,1)))
         else
            do ifile=1,nfiles
               read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3,v4
               res(ifile-1) = v3
            enddo

            call getpdfuncertainty(res(0:nmempdf_end),central,errplus
     $           ,errminus,errsymm)

            write(51,'(4(1x,d14.8))') v1,v2,central,errsymm
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
