c     vegas with double precision
c...............
      module integration
      
      logical readin,writeout
      character*72 ingridfile,outgridfile, outgridtopfile
      common/gridinfo_logic/readin,writeout
      common/gridinfo_char/ingridfile,outgridfile,outgridtopfile

      contains

      SUBROUTINE vegas(region,ndim,fxn,init,ncall,itmx,nprn,
     >     tgral,sd,chi2a)
      INTEGER init,itmx,ncall,ndim,nprn,NDMX,MXDIM
      REAL*8 tgral,chi2a,sd,region(2*ndim),fxn,ALPH,TINY
      PARAMETER (ALPH=1.5,NDMX=50,MXDIM=11,TINY=1.d-30)
      EXTERNAL fxn
C     
      INTEGER i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(MXDIM),kg(MXDIM)
      double precision
     1     calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,
     >     d(NDMX,MXDIM),di(NDMX,MXDIM),dt(MXDIM),dx(MXDIM),
     >     r(NDMX),x(MXDIM),xi(NDMX,MXDIM),xin(NDMX),xval
      integer vegas_ncall
      common/vegas_ncall/ vegas_ncall
      DOUBLE PRECISION schi,si,swgt
      COMMON /ranno/ idum
      common/vegas_wgt/wgt,swgt
      integer ilast
      common/last_integ/ilast
      logical set_ilast
      SAVE
      if(init.le.0)then
         mds=1
         ndo=1
         do 11 j=1,ndim
            xi(1,j)=1
 11      enddo
      endif
      if (init.le.1) then
         if(ilast.ne.1) goto 29
         
 29      si=0
         swgt=0
         schi=0
      endif
      if (init.le.2)then
         nd=NDMX
         ng=1
         if(mds.ne.0)then
            ng=(ncall/2.d0+0.25d0)**(1.d0/ndim)
            mds=1
            if((2*ng-NDMX).ge.0)then
               mds=-1
               npg=ng/NDMX+1
               nd=ng/npg
               ng=npg*nd
            endif
         endif
         k=ng**ndim
         npg=max(ncall/k,2)
         calls=npg*k
         vegas_ncall=calls
         dxg=1.d0/ng
         dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.d0)
         xnd=nd
         dxg=dxg*xnd
         xjac=1.d0/calls
         do 12 j=1, ndim
            dx(j)=region(j+ndim)-region(j)
            xjac=xjac*dx(j)
 12      enddo

c---  read-in grid if necessary
         if (readin) then
            print*, ingridfile, '*'
            open(unit=13,file=ingridfile,status='unknown')
            write(6,*)'****************************************************'
            write(6,*)'* Reading in vegas grid from ',ingridfile,' *'
            write(6,*)'****************************************************'
            call flush(6)
            do j=1,ndim
               read(13,203) jj,(xi(i,j),i=1,nd)
            enddo
            close(13)
            ndo=nd
            readin=.false.
         endif


         if(nd.ne.ndo)then
            do 13 i=1,nd
               r(i)=1.d0
 13         enddo
            do 14 j=1,ndim
               call rebin(ndo/xnd,nd,r,xin,xi(1,j))
 14         enddo
            ndo=nd
         endif
         if(nprn.ge.0) write(*,200) ndim,calls,it,itmx,nprn,
     *        ALPH,mds,nd,(j,region(j),j,region(j+ndim),j=1,ndim)
      endif
      do 28 it=1,itmx
         if (ilast.eq.1.and.it.lt.itmx) then
            ilast=0
            set_ilast=.true.
         endif
         if (it.eq.itmx.and.set_ilast) then
            ilast = 1
         endif
         ti=0.d0
         if (ilast.ne.1.or.it.ne.itmx) goto 30
 30      tsi=0.d0
         do 16 j=1,ndim
            kg(j)=1
            do 15 i=1,nd
               d(i,j)=0.d0
               di(i,j)=0.d0
 15         enddo
 16      enddo
 10      continue
         fb=0.d0
         if (ilast.ne.1.or.it.ne.itmx) goto 31
 31      f2b=0.d0
         do 19 k=1,npg
            wgt=xjac
            do 17 j=1,ndim
               xn=(kg(j)-ran2(idum))*dxg+1.d0
               ia(j)=max(min(int(xn),NDMX),1)
               if(ia(j).gt.1)then
                  xo=xi(ia(j),j)-xi(ia(j)-1,j)
                  rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
               else
                  xo=xi(ia(j),j)
                  rc=(xn-ia(j))*xo
               endif
               x(j)=region(j)+rc*dx(j)
               wgt=wgt*xo*xnd
 17         enddo
            f=wgt*fxn(x,wgt)
            f2=f*f
            fb=fb+f

            if (ilast.ne.1.or.it.ne.itmx) goto 32


 32         f2b=f2b+f2
            do 18 j=1,ndim
               di(ia(j),j)=di(ia(j),j)+f
               if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
 18         enddo
 19      enddo
         f2b=sqrt(f2b*npg)
         f2b=(f2b-fb)*(f2b+fb)
         if (f2b.le.0) f2b=TINY
         ti=ti+fb
         if (ilast.ne.1.or.it.ne.itmx) goto 33
 33      tsi=tsi+f2b
         if(mds.lt.0)then
            do 21 j=1,ndim
               d(ia(j),j)=d(ia(j),j)+f2b
 21         enddo
         endif
         do 22 k=ndim,1,-1
            kg(k)=mod(kg(k),ng)+1
            if(kg(k).ne.1) goto 10
 22      enddo
         tsi=tsi*dv2g
         wgt=1.d0/tsi

         si=si+dble(wgt)*dble(ti)
         if (ilast.ne.1.or.it.ne.itmx) goto 34
 34      schi=schi+dble(wgt)*dble(ti)**2
         swgt=swgt+dble(wgt)
         tgral=si/swgt
         if (ilast.ne.1.or.it.ne.itmx) goto 35
 35      chi2a=max((schi-si*tgral)/(it-.99d0),0.d0)
         sd=sqrt(1./swgt)
         tsi=sqrt(tsi)
         if(nprn.ge.0)then
            write(*,201) it,ti,tsi,tgral,sd,chi2a
            if(nprn.ne.0)then
               do 23 j=1,ndim
                  write(*,202) j,(xi(i,j),di(i,j),
     *                 i=1+nprn/2,nd,nprn)
 23            enddo
            endif
         endif
         do 25 j=1,ndim
            xo=d(1,j)
            xn=d(2,j)
            d(1,j)=(xo+xn)/2.d0
            dt(j)=d(1,j)
            do 24 i=2,nd-1
               rc=xo+xn
               xo=xn
               xn=d(i+1,j)
               d(i,j)=(rc+xn)/3.d0
               dt(j)=dt(j)+d(i,j)
 24         enddo
            d(nd,j)=(xo+xn)/2.d0
            dt(j)=dt(j)+d(nd,j)
 25      enddo
         do 27 j=1,ndim
            rc=0.d0
            do 26 i=1,nd
               if(d(i,j).lt.TINY) d(i,j)=TINY
               r(i)=((1.d0-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
               rc=rc+r(i)
 26         enddo
            call rebin(rc/xnd,nd,r,xin,xi(1,j))
 27      enddo
 28   enddo

c---  write-out grid if necessary
      if (writeout) then
         open(unit=11,file=outgridfile,status='unknown')
         open(unit=111,file=outgridtopfile,status='unknown')
         write(6,*)'****************************************************'
         write(6,*)'* Writing out vegas grid to ',outgridfile,'  *'
         write(6,*)'****************************************************'
         call flush(6)
         do j=1,ndim
            write(11,203) jj,(xi(i,j),i=1,nd)
            write(111,*) '# XDIM index  ', j-1
            do i=1,nd
               xval=(i*1d0)/nd
               write(111,'(D16.8,D16.8)') xval, xi(i,j)
            enddo
            write(111,*) ''
            write(111,*) ''
         enddo
         close(11)
         close(111)
      endif

      return
 200  FORMAT(/' input parameters for vegas: ndim=',i3,' ncall=',f8.0
     *     /28x,' it=',i5,' itmx=',i5
     *     /28x,' nprn=',i3,' alph=',f5.2/28x,' mds=',i3,' nd=',i4
     *     /(30x,'xl(',i2,')= ',g11.4,' xu(',i2,')= ',g11.4))
 201  FORMAT(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2
     *     /' all iterations: integral =',g14.7,'+/- ',g9.2,
     *     ' chi**2/it''n =',g9.2)
 202  FORMAT(/' data for axis ',I2/' X delta i',
     *     ' x delta i ',' xdelta i ',
     *     /(1x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
 203  FORMAT(/(5z16))
      end subroutine

c...............................................................
      SUBROUTINE rebin(rc,nd,r,xin,xi)
      INTEGER nd
      REAL*8 rc,r(*),xi(*),xin(*)
      INTEGER i,k
      REAL*8 dr,xn,xo
      k=0
      xn=0.d0
      dr=0.d0
      do 11 i=1,nd-1
 1       if(rc.gt.dr)then
            k=k+1
            dr=dr+r(k)
            xo=xn
            xn=xi(k)
            goto 1
         endif
         dr=dr-rc
         xin(i)=xn-(xn-xo)*dr/r(k)
 11   enddo
      do 12 i=1,nd-1
         xi(i)=xin(i)
 12   enddo
      xi(nd)=1.d0
      return
      end subroutine
c...............................................................
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      enddo
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function
c....................................................................

      end module integration
