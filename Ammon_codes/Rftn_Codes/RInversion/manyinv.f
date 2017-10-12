C Output from Public domain Ratfor, version 1.0
      program manyinv
      parameter(nlmax = 45, ntmax = 520, nsmax = 2, ndat = ntmax*nsmax+2
     **nlmax)
      dimension alpha(nlmax),beta(nlmax),rho(nlmax),thiki(nlmax),h(nlmax
     *)
      dimension pert(nlmax),alphai(nlmax),betai(nlmax),rhoi(nlmax)
      character*32 modela,title
      integer inunit,ounit,ifile
      logical porsv(nsmax)
      common /seismo/ seis(ntmax,nsmax), dt(nsmax), dura(nsmax), dly(nsm
     *ax),gauss(nsmax),p(nsmax),nt(nsmax),porsv
      common /imodel/alpha,beta,thiki,rho,nlyrs
      common /innout/ inunit,ounit
      character*24 todays_date
      real tfraction
      real fmin
      integer npasses
      logical hpfilter, yesno
      common /filter/ fmin, npasses, hpfilter
      inunit = 5
      ounit = 6
      ifile = 8
      seed = 0.5
      fmin = 0.03
      npasses = 2
      hpfilter = .false.
      call rand0(seed)
      write(ounit,'(/)')
      write(ounit,*) '**************************************************
     *********'
      write(ounit,*)'manyinv - Receiver function inversion program.'
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
      pert(i) = 0.0
      alphai(i) = 0.0
      alpha(i) = 0.0
23000 continue
23001 continue
      do23002 i = 1,nsmax
      porsv(i) = .true.
23002 continue
23003 continue
      write(ounit,*) ' '
      write(ounit,*)'input velocity model:'
      read(inunit,'(a)')modela
      write(ounit,*)'maximum perturbation in km/sec'
      read(inunit,*) pertmax
      write(ounit,*)'Velocity to cut perturbing off'
      read(inunit,*) vcut
      write(ounit,*)'Maximum perturbation for random component'
      write(ounit,*)'in percent of the maximum perturbation input'
      write(ounit,*)'above (10. -> 10%)'
      read(inunit,*) rpercent
      write(ounit,*) ' '
      write(ounit,*)'Enter the max number of iterations per inversion'
      read(inunit,*) maxiter
      write(ounit,*)'Enter the smoothing trade-off parameter'
      read(inunit,*) sigjmpb
      write(ounit,*) ' '
      write(ounit,*)'Initial models are generated in 2 loops,'
      write(ounit,*)'when you enter the number of inversions'
      write(ounit,*)'you will actually get 4 times that number'
      write(ounit,*)'adjust your choice to compensate for this.'
      write(ounit,*) ' '
      write(ounit,*) 'So how many inversions do you want me to do?'
      read(inunit,*) numinv
      write(ounit,*) ' '
      write(ounit,*)'Enter Singular Value truncation fraction'
      read(inunit,*) tfraction
      hpfilter = yesno('Apply a high-pass filter to waveforms? ')
      if(.not.(hpfilter))goto 23004
      write(ounit,*)'Enter the corner frequency.'
      read(inunit,*) fmin
      write(ounit,*) 'Enter the number of filter passes (1 or 2).'
      read(inunit,*) npasses
23004 continue
      rpert = rpercent/100.0
      call getseis(ns,seis,ntmax,nsmax,dt,dura,dly,gauss,p,nt,porsv)
      open(unit=ifile,file=modela)
      rewind=ifile
      read(ifile,100)nlyrs,title
100   format(i3,1x,a32)
      do23006 i1 = 1,nlyrs 
      read(ifile,110)idum,alphai(i1),betai(i1),rhoi(i1),thiki(i1),dum1,d
     *um2,dum3,dum4,dum5
23006 continue
23007 continue
110   format(i3,1x,9f8.4)
      close(unit=ifile)
      tdpth = 0.
      nlc = nlyrs
      iflag = 0
      do23008 i2 = 2,nlyrs 
      itemp = i2 - 1
      tdpth = tdpth + thiki(itemp)
      h(i2) = tdpth
      if(.not.(alphai(i2).le.vcut))goto 23010
      cthick = tdpth+thiki(i2+1)
      iflag = 1
      nlc = i2
23010 continue
23008 continue
23009 continue
      if(.not.(iflag .eq. 0))goto 23012
      cthick = thiki(nlyrs)
23012 continue
      h(1) = 0.
      r1 = -1.
      do23014 iouter = 1, 4 
      r2 = 0.0
      r3 = 1.0
      do23016 inner = 1, numinv 
      r2 = r2 + float(inner-1)/float(numinv)
      a2 = -(r1 + r2 + r3)
      a1 = r1*r2 + r1*r3 + r2*r3
      a0 = -(r1 * r2 * r3)
      amax = 0.0
      do23018 i5 = 1,nlc 
      z = h(i5)/cthick
      pert(i5) = cubic(z,a2,a1,a0)
      if(.not.(amax .le. abs(pert(i5))))goto 23020
      amax = abs(pert(i5))
23020 continue
23018 continue
23019 continue
      anorm = pertmax/amax
      do23022 i6 = 1,nlyrs 
      call rand0(seed)
      randpart = 2.*(seed - 0.5) * pertmax * rpert
      alpha(i6) = alphai(i6) + pert(i6) * anorm + randpart
      beta(i6) = alpha(i6)/1.732050808
      rho(i6) = 0.32 * alpha(i6) + 0.77
23022 continue
23023 continue
      invnum = float(iouter-1)*numinv + inner
      write(ounit,*) '==================================================
     *========'
      write(ounit,*) '==================================================
     *========'
      write(ounit,*) '    Inversion Number: ',invnum
      write(ounit,*) '==================================================
     *========'
      write(ounit,*) '==================================================
     *========'
      call jinv(sigjmpb,maxiter,ns,invnum,tfraction)
23016 continue
23017 continue
      r1 = r1 +.5*float(iouter)
23014 continue
23015 continue
      stop
      end
      subroutine rand0(x)
      data k,j,m,rm/5701,3612,566927,566927.0/
      ix=int(x*rm)
      irand=mod(j*ix+k,m)
      x=(real(irand)+.5)/rm
      return
      end
      real function cubic(z,a2,a1,a0)
      cubic = a0+z*(a1+z*(a2+z))
      end
