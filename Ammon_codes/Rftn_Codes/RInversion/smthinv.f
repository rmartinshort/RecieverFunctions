C Output from Public domain Ratfor, version 1.0
      program smthinv
      parameter(nlmax = 45, ntmax = 520, nsmax = 2, ndat = ntmax*nsmax+2
     **nlmax)
      dimension alpha(nlmax),beta(nlmax),rho(nlmax),thiki(nlmax)
      character*32 modela,title
      real minsigma,maxsigma,dsigma
      integer inunit,ounit,oun2,icount
      logical porsv(nsmax)
      character*24 todays_date
      common /seismo/ seis(ntmax,nsmax), dt(nsmax), dura(nsmax), dly(nsm
     *ax),gauss(nsmax),p(nsmax),nt(nsmax),porsv
      common /imodel/alpha,beta,thiki,rho,nlyrs
      common /innout/ inunit,ounit
      real tfraction
      real fmin
      integer npasses
      logical hpfilter, yesno
      common /filter/ fmin, npasses, hpfilter
      inunit = 5
      ounit = 6
      oun2 = 8
      write(ounit,'(/)')
      write(ounit,*) '**************************************************
     *********'
      write(ounit,*)'smthinv - Receiver function inversion program.'
      write(ounit,*)'          VERSION 2.1 July 1997'
      write(ounit,*)'    Charles J. Ammon and George Randall.'
      write(ounit,*)'Additional routines by George Zandt and Tom Owens.'
      write(ounit,*) '**************************************************
     *********'
      call fdate(todays_date)
      write(ounit,*) 'Inversion run on: ',todays_date
      write(ounit,*) '**************************************************
     *********'
      write(ounit,*)'Maximum Number of points in each waveform = 512'
      write(ounit,*) '**************************************************
     *********'
      do23000 i = 1,nlmax 
      alpha(i) = 0.
      beta(i) = 0.
      rho(i) = 0.
      thiki(i) = 0.
23000 continue
23001 continue
      do23002 i = 1,nsmax
      porsv(i) = .true.
23002 continue
23003 continue
      write(ounit,*)'input velocity model:'
      read(inunit,'(a)')modela
      write(ounit,*)'Enter the max number of iterations per inversion'
      read(inunit,*) maxiter
      write(ounit,*)'Enter the minimum smoothing trade-off parameter'
      read(inunit,*) minsigma
      write(ounit,*)'Enter the maximum smoothing trade-off parameter'
      read(inunit,*) maxsigma
      write(ounit,*)'Enter Singular Value truncation fraction'
      read(inunit,*) tfraction
      hpfilter = yesno('Apply a high-pass filter to waveforms? ')
      if(.not.(hpfilter))goto 23004
      write(ounit,*)'Enter the corner frequency.'
      read(inunit,*) fmin
      write(ounit,*) 'Enter the number of filter passes (1 or 2).'
      read(inunit,*) npasses
23004 continue
      call getseis(ns,seis,ntmax,nsmax,dt,dura,dly,gauss,p,nt,porsv)
      icount = 1
      dsigma = (maxsigma - minsigma)/10
23006 if(.not.(sigjmpb .le. maxsigma))goto 23007
      sigjmpb = minsigma + (icount-1)*dsigma
      write(ounit,*) ' '
      write(ounit,*) '**************************************************
     *********'
      write(ounit,*)'Smoothness trade-off parameter = ', sigjmpb
      write(ounit,*) '**************************************************
     *********'
      write(ounit,*) ' '
      open(unit=oun2,file=modela)
      rewind=oun2
      read(oun2,100)nlyrs,title
100   format(i3,1x,a32)
      do23008 i1 = 1,nlyrs 
      read(oun2,110)idum,alpha(i1),beta(i1),rho(i1),thiki(i1),dum1,dum2,
     *dum3,dum4,dum5
23008 continue
23009 continue
110   format(i3,1x,9f8.4)
      close(unit=oun2)
      invnum = icount
      call jinv(sigjmpb,maxiter,ns,invnum,tfraction)
      icount = icount + 1
      goto 23006
23007 continue
      stop
      end
