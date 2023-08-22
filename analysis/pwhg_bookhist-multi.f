      subroutine setupmulti(n)
      implicit none
      integer n
      include 'pwhg_bookhist-multi.h'
      real * 8 weirdnum
      common/c_setupmulti/weirdnum
      save /c_setupmulti/
c if this is set we are sure that setupmulti was called
      weirdnum=317d0/12345d0
      if(n.gt.maxmulti) then
         write(*,*) ' ************** ERROR **************'
         write(*,*) ' multi-weight histogramming requested with',n,
     1 'weights > ',maxmulti
         write(*,*) ' increase maxmulti in pwhg-bookhist-multi.h'
         call exit(1)
      endif
      nmulti=n
      end

      subroutine bookupeqbins(string,binsize,xlow,xhigh)
      implicit none
      character *(*) string
      real * 8 binsize,xlow,xhigh
      include 'pwhg_bookhist-multi.h'
      real * 8 xx(maxbins+1)
      integer k
      xx(1)=xlow
      do k=2,maxbins+1
         xx(k)=xx(k-1)+binsize
c avoid funny bin edges near zero ...
         if(abs(xx(k)/binsize).lt.1d-6) xx(k) = 0
         if(xx(k)-(xhigh-binsize/1e4).gt.0) goto 10
      enddo
      write(*,*) 'bookupeqbins: too many bins in hist ',string
      call exit(-1)
 10   continue
      if((xx(k)-xhigh)/binsize.gt.1e-4) then
         write(*,*) 'upper limit incompatible with bin size'
         write(*,*) 'replacing ',xhigh,' with ',xx(k)
         write(*,*) ' in histogram ',string
      endif
      call bookup(string,k-1,xx)
      end


      subroutine bookup(string,n,x)
c Books up a histogram characterized by the tag string <string>,
c with n bins. The array x(n+1) are the bins endpoints,
c x(i) is the low extreme of bin i.
      implicit none
      character *(*) string
      integer n
      real * 8 x(n+1)
      include 'pwhg_bookhist-multi.h'
      integer j,k
      integer indexhist
      real * 8 weirdnum
      common/c_setupmulti/weirdnum
      save /c_setupmulti/
      logical ini
      data ini/.true./
      save ini

c We assume that this routine is always called first when the package is used
      if(ini) then
         if(.not.weirdnum.eq.317d0/12345d0) then
c setupmulti was not called! setup default value!
            call setupmulti(1)
         endif
         ini=.false.
      endif

      if(n.gt.maxbins) then
         write(*,*) ' maximum number of bins=',maxbins
         write(*,*) ' asked for ',n
         call exit(-1)
      endif
c indexhist(string) returns the histogram index if a histogram
c with tag string was already booked, otherwise it books a new histogram,
c and returns minus the value of its index
      j=-indexhist(string)
      if(j.lt.0) then
         write(*,*) 'Histogram ',string,' already booked'
         call exit(-1)
      endif
      do k=1,n+1
         xhistarr(k,j)=x(k)
      enddo
c y and err values go from 0 to n+1, 0 being the underflow and n+1
c the overflow.
      do k=0,n+1
         yhistarr(:,k,j)=0
         yhistarr1(:,k,j)=0
         errhistarr1(:,k,j)=0
         yhistarr2(:,k,j)=0
         errhistarr2(:,k,j)=0
         nhits(k,j)=0
      enddo
      nbins(j)=n
      ient1(j)=0
      end
      

      function indexhist(string)
      implicit none
      character * (*) string
      include 'pwhg_bookhist-multi.h'
      integer indexhist
      integer j,khist
      if(string.eq.' ') then
         write(*,*) ' indexhist: error, empty name'
         call exit(-1)
      endif
      khist = -1
      do j=1,jhist
         if(stringhist(j).eq.' ') then
            khist = j
         endif
         if(stringhist(j).eq.string) then
            indexhist=j
            return
         endif
      enddo
      if(khist.gt.0) then
         goto 999
      endif
      if(jhist.eq.nhist) then
         write(*,*) ' no more rooms for histograms'
         write(*,*) ' Histogram "',string,'" cannot be booked'
         call exit(-1)
      endif
      jhist=jhist+1
      khist = jhist
 999  stringhist(khist)=trim(adjustl(string))
      if(stringhist(khist).ne.trim(adjustl(string))) then
         write(*,*) ' Histogram string "',string,'" too long'
         call exit(-1)
      endif
c the negative sign indicates a new histogram
      indexhist=-khist
      end


      subroutine deletehist(string,iret)
      implicit none
      character * (*) string
      integer iret
      include 'pwhg_bookhist-multi.h'
      integer j
      if(string.eq.' ') then
         write(*,*) ' deletehist: error, empty name'
         call exit(-1)
      endif
      do j=1,jhist
         if(stringhist(j).eq.string) then
            stringhist(j)=' '
            return
            iret = 0
         endif
      enddo
      write(*,*) ' deletehist: histogram '//string
     1     //' not present'
      iret = -1
      end

      subroutine filld(string,xval,weight)
      implicit none
      character *(*) string
      include 'pwhg_bookhist-multi.h'
      real * 8 xval,weight(1:nmulti)
      integer j,k,indexhist
      j=indexhist(string)
      if(j.lt.0) then
         write(*,*) ' histogram "',string,'" was not booked'
         call exit(-1)
      endif
c underflow
      if(xval.lt.xhistarr(1,j)) then
         yhistarr(1:nmulti,0,j)=yhistarr(1:nmulti,0,j)+weight
         nhits(0,j)=nhits(0,j)+1
         return
      else
         do k=1,nbins(j)
            if(xval.lt.xhistarr(k+1,j)) then
               yhistarr(1:nmulti,k,j)=yhistarr(1:nmulti,k,j)+weight/
     1              (xhistarr(k+1,j)-xhistarr(k,j))
               nhits(k,j)=nhits(k,j)+1
               return
            endif
         enddo
      endif
c overflow
      yhistarr(1:nmulti,nbins(j)+1,j)=yhistarr(1:nmulti,nbins(j)+1,j)
     1     +weight
      end


      subroutine inihists
      implicit none
      include 'pwhg_bookhist-multi.h'
      jhist=0
      end

      subroutine pwhgtopout(filename)
      use parameters
      implicit none
      character * (*) filename
      include 'pwhg_bookhist-multi.h'
      integer k,j,iun,l
      character * 3 cl
      call newunit(iun)
      do l=1,nmulti
         if(nmulti.eq.1) then
            if(outname.ne."") filename = outname
            open(unit=iun,file=trim(adjustl(filename))//'.top',
     1           status='unknown')
         else
            write(cl,'(i3)') l
!            open(unit=iun,file=trim(adjustl(filename))//'-W'//
!     1           trim(adjustl(cl))//'.top',status='unknown')
            open(unit=iun,file=trim(adjustl(filename))
     $           //trim(adjustl(scalestr(l)))//'.top',status='unknown')
         endif
         do j=1,jhist
            write(iun,'(a,i3)')'# '//trim(adjustl(stringhist(j)))//
     1           ' index ',j-1
            do k=1,nbins(j)
               write(iun,'(4(1x,e14.8))') xhistarr(k,j),
     1              xhistarr(k+1,j),yhistarr2(l,k,j),errhistarr2(l,k,j)
            enddo
            write(iun,*)
            write(iun,*)
         enddo
         close(iun)
      enddo
      end


      subroutine pwhgaccumup
c values histogrammed so far are transferred to array yhistarr1,
c and the square of the values are transferred to array errhistarr1.
c yhistarr is zeroed. The index ient1 is increased by one unit.
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer j,k
      do j=1,jhist
         do k=0,nbins(j)+1
            yhistarr1(1:nmulti,k,j)=yhistarr1(1:nmulti,k,j)
     1           +yhistarr(1:nmulti,k,j)
            errhistarr1(1:nmulti,k,j)=errhistarr1(1:nmulti,k,j)
     1           +yhistarr(1:nmulti,k,j)**2
            yhistarr(1:nmulti,k,j)=0
         enddo
         ient1(j)=ient1(j)+1
      enddo
      end

      subroutine pwhgsetout
c provides a snapshot of the current result of the
c analysis, leaving the yhistarr1 and errhistarr1 unchanged.
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer j,k
      real *8 xxx,sum(1:nmulti),sumsq(1:nmulti)
      do j=1,jhist
         xxx=1d0/ient1(j)
         do k=0,nbins(j)+1
            sum=yhistarr1(1:nmulti,k,j)
            sumsq=errhistarr1(1:nmulti,k,j)
            yhistarr2(1:nmulti,k,j)=xxx*sum
            errhistarr2(1:nmulti,k,j)=sqrt(xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      end

      subroutine pwhgaddout
c accumulates the results obtained so far in yhistarr2 and errhistarr2.
c It zeroes yhistarr1 and errhistarr1. To be used if we compute
c a cross section with several contributions.
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer j,k
      real *8 xxx,sum(1:nmulti),sumsq(1:nmulti)
      do j=1,jhist
         xxx=1d0/ient1(j)
         do k=0,nbins(j)+1
            sum=yhistarr1(1:nmulti,k,j)
            sumsq=errhistarr1(1:nmulti,k,j)
            yhistarr2(1:nmulti,k,j)=yhistarr2(1:nmulti,k,j)+xxx*sum
            errhistarr2(1:nmulti,k,j)=sqrt(errhistarr2(1:nmulti,k,j)**2+
     1           xxx**2*abs(sumsq-sum**2*xxx))
         enddo
      enddo
      do j=1,jhist
         do k=0,nbins(j)+1
            yhistarr1(1:nmulti,k,j)=0
            errhistarr1(1:nmulti,k,j)=0
         enddo
         ient1(j)=0
      enddo
      end

      subroutine resethists
c Set the yhistarr* and errhistarr* arrays to zero for use when we do
c multiple runs within one call (ie for scale and pdf variations)
      implicit none
      include 'pwhg_bookhist-multi.h'
      integer j,k
      real *8 xxx,sum(1:nmulti),sumsq(1:nmulti)

      yhistarr    = 0d0
      yhistarr1   = 0d0
      yhistarr2   = 0d0
      errhistarr1 = 0d0
      errhistarr2 = 0d0
      ient1       = 0


      end


      subroutine newunit(iun)
      implicit none
      integer iun
      logical ok
      integer j
c units 97 and 99 are used for lhe files;
c keep reserved.
      do j=10,96
         inquire(unit=j,opened=ok)
         if(.not.ok) then
            iun=j
            return
         endif
      enddo
      end
