c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine define_histograms
         use mod_parameters
         implicit none
         include 'pwhg_bookhist-multi.h'
         integer xnbins, Q2nbins
         
         call bookupeqbins('sigincl',1d0,0d0,1d0)
         call bookupeqbins('E/Q/2macro',10d0,0d0,sqrts)
         call bookupeqbins('logE/Q/2macro',0.2d0,-5.5d0,4.5d0)
         call bookupeqbins('etamacro',0.5d0,-40d0,15d0)
          
      end
      
      subroutine user_analysis(n,dsig,x,y,Q2)
         use mod_parameters
         use mod_analysis
         implicit none
         real * 8 dsig(maxscales)

         integer n
         double precision x, y, Q2
         integer maxjets
         parameter (maxjets=3)
         double precision ppartons(0:3,maxjets),pj(0:3,maxjets),q(0:3)
         double precision ptj1, etaj1
         real * 8 sqrtQ2,sbeams,kT2, palg, kt, eta, Ej, etaj

         integer i,npartons,njets,beam_sign, imacrojet, macro_jet_index

         double precision Ecut, Ecur
         logical Ecur_lt_Ecut, CutDIS_Ecur

         call filld('sigincl',0.5d0,dsig)

         sqrtQ2 = sqrt(Q2)
         Ecut = 0d0! sqrtQ2/10d0

         npartons = n - 3
c        Fill the parton momenta array with breit frame partons
         q = pbreit(:,1) - pbreit(:,3)
         ppartons(:,1:npartons) = pbreit(:,4:n)

         Ecur_lt_Ecut = CutDIS_Ecur(npartons, ppartons(:,:), q, Ecut, Ecur)

         if(Ecur_lt_Ecut) return

         palg = 1d0 ! 0: C/A, 1: kT like
         kT2 = 900d0 ! in GeV^2
         beam_sign = sign(1d0,pbreit(2,1)) ! proton beam direction
         !print*, 'Clustering jets with DIS Cambridge algorithm, kT2=', kT2, ' palg=', palg, ' beam_sign=', beam_sign
         call buildjets(npartons,ppartons,kT2,palg,beam_sign,pj,njets)

         !print*, '--- Event analysis ---' 
         !print*, 'Q, x, y:', sqrtQ2, x, y
         !print*, 'Number of jets found:', njets
         !do i =1,njets
         !   print*, 'Jet', i, 'E, pT, eta:', pj(0,i), kt(pj(:,i)), eta(pj(:,i))
         !enddo
         !print*, ''

         ! Select the macro jet as the jet with the largest z component
         imacrojet = macro_jet_index(pj,njets,beam_sign)
         etaj = eta(pj(:,imacrojet))
         Ej = pj(0,imacrojet)

         call filld('E/Q/2macro',Ej/(sqrtQ2/2d0),dsig)
         call filld('logE/Q/2macro',log10(Ej/(sqrtQ2/2d0)),dsig)
         call filld('etamacro',max(-5d100+1d-10,etaj),dsig)
      end

c     For Ecut < Q*(sqrt(2)-1) this selects only events with no
c     particles in the current hemisphere For Ecut < 0.5*Q, only real
c     radiation is involved, hence the cut selects DIS two-jet events
      logical function CutDIS_Ecur(n,p,q,Ecut,Ecur) 
      implicit none
      integer n
      double precision p(0:3,n),q(0:3), Ecut
!----------------------------------------
      integer :: i
      double precision :: Ecur
      
      Ecur = 0d0
      do i=1,n
          if (dot_product(p(1:3,i),q(1:3))>0d0) then
            Ecur = Ecur + p(0,i)
         end if
      end do
      cutDIS_Ecur = Ecur<Ecut    
      
      end function CutDIS_Ecur

      subroutine buildjets(n, pin, kT2, palg, beam_sign, pj, njets)
         implicit none
         integer n
         double precision pin(0:3,n), kT2, palg
         integer maxtrack,maxjet
         parameter (maxtrack=3,maxjet=3)

   !     Output
         double precision pj(0:3,maxjet)
   !     Internal
         integer mu, njets, ntracks, ijet, j, jetvec(maxtrack)
         double precision pjet(4,maxjet), ptmin
         double precision ptrack(4,maxtrack)
         integer beam_sign

         ptrack = 0d0
         pjet = 0d0
         njets=0
         ntracks = n
         pj = 0d0
         ptmin = 0d0
         jetvec = 0

   !     Fast jet needs 0 indexed tracks
         ptrack(4,1:n)=pin(0,1:n)
         do mu=1,3
            ptrack(mu,1:n)=pin(mu,1:n)
         enddo

         call fastjetdiscambridge(ptrack,ntracks,kT2,palg,beam_sign,pjet,njets)

         if(njets.gt.3.or.njets.lt.1) then
            print*, 'njets out of bounds!!', njets
            stop
         endif

   !     Back to 1:4 index
         do ijet=1,njets
            do mu=1,3
               pj(mu,ijet)=pjet(mu,ijet)
            enddo
            pj(0,ijet)=pjet(4,ijet)
         enddo
      end

      double precision function kt(p)
         implicit none
         double precision p(0:3)
         
         kt = sqrt(max(p(1)**2 + p(2)**2,0d0))
      end

      function eta(p)
         implicit none
         real * 8 eta, p(0:3), normp, norm

         normp = norm(p)
         if(normp.gt.p(3)) then
            eta = 0.5d0 * log((normp + p(3)) / (normp - p(3)))
         else
            eta = sign(1d100,p(3)) 
         endif
      end

      function norm(p)
         implicit none
         real * 8 norm, p(0:3)

         norm = p(1)**2 + p(2)**2 + p(3)**2
         norm = sqrt(max(0d0,norm))
      end

      function getrapidity0(p)
         implicit none
         real * 8 p(0:3),getrapidity0
         getrapidity0=0.5d0*log((p(0)+p(3))/(p(0)-p(3)))
      end

      subroutine getrapidity(p,y)
         implicit none
         real * 8 p(4),y
         y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      end

      function azi(p)
         implicit none
         real * 8 pi
         parameter(pi = 3.141592653589793D0)
         real * 8 azi,p(0:3)
         azi = atan(p(2)/p(1))
         if (p(1).lt.0d0) then
            if (azi.gt.0d0) then               
               azi = azi - pi
            else
               azi = azi + pi
            endif
         endif    
      end

      integer function macro_jet_index(pj,njets,beam_sign)
         implicit none
         integer maxjets
         parameter (maxjets = 3)
         double precision pj(0:3,maxjets)
         integer njets, beam_sign
         integer i, imacro
         double precision longmax, long
         longmax = -1d100
         do i=1,njets
            long = pj(0,i) - pj(3,i)*beam_sign
            if (long.gt.longmax) then
               longmax = long
               imacro = i
            endif
         enddo
         macro_jet_index = imacro
      end
