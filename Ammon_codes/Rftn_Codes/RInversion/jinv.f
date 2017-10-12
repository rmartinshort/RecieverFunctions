C Output from Public domain Ratfor, version 1.0
      subroutine jinv(sigjmpb,maxiter,ns,invnum,tfraction)
      parameter ( nlmax = 45, ntmax = 520, nsmax = 2, ndat = ntmax*nsmax
     *+2*nlmax)
      real a(ndat,2*nlmax), b(ndat)
      real s(3*nlmax),wk(4*nlmax)
      logical porsv(nsmax)
      integer nt(nsmax)
      common /seismo/ seis(ntmax,nsmax), dt(nsmax), dura(nsmax), dly(nsm
     *ax), gauss(nsmax), p(nsmax),nt,porsv
      real fmin
      integer npasses
      logical hpfilter
      common /filter/ fmin, npasses, hpfilter
      real aa(ntmax,nlmax),rms(nsmax)
      real syn(ntmax,nsmax)
      real sol(2*nlmax)
      real tfraction, sig_power(nsmax),misfit(nsmax)
      real perta(nlmax), pertb(nlmax), pertr(nlmax)
      real vp_over_vs(nlmax),pratio(nlmax)
      logical pon(nlmax,6)
      common /imodel/ alpha(nlmax), beta(nlmax), thk(nlmax), rho(nlmax),
     *nlyrs
      logical happy, jumping, yesno
      integer inunit,ounit,iter,maxiter
      common /innout/ inunit,ounit
      inunit=5
      ounit=6
      iter=0
      jumping = .true.
      sigjmpa = 0.
      do23000 i = 1, nlyrs-1 
      pon(i,1) = .true. 
      pon(i,2) = .true. 
      pon(i,3) = .true.
23000 continue
23001 continue
      pon(nlyrs,1) = .true. 
      pon(nlyrs,2) = .true. 
      pon(nlyrs,3) = .true.
      happy = .false.
      write(*,*) 'Initial Model Vp over Vs Ratio'
      write(*,*) 'Layer   Vp/Vs       Poissons Ratio'
      do23002 ilyr=1,nlyrs
      vp_over_vs(ilyr) = alpha(ilyr)/beta(ilyr)
      pratio(ilyr) = vpovs_to_pr(vp_over_vs(ilyr))
      write(*,'(i5,2x,f10.7,8x,f5.3)') ilyr, vp_over_vs(ilyr),pratio(ily
     *r)
23002 continue
23003 continue
      if(.not.(hpfilter))goto 23004
      do23006 iseis = 1, ns 
      call hpbufilter(seis(1,iseis),nt(iseis),dt(iseis),fmin,npasses)
23006 continue
23007 continue
23004 continue
      call putseis(ns, seis, ntmax, nsmax, dt, dura, dly, gauss, p, nt, 
     *porsv, iter, invnum)
      do23008 iseis = 1, ns 
      sig_power(iseis) = 0.0
      do23010 j = 1, nt(iseis)
      sig_power(iseis) = sig_power(iseis) + seis(j,iseis)*seis(j,iseis)
23010 continue
23011 continue
      if(.not.(sig_power(iseis) .eq. 0.0))goto 23012
      write(*,*) 'Signal ',iseis,' has zero power, stopping execution.'
      stop
23012 continue
23008 continue
23009 continue
      do23014 i = 1, ndat 
      do23016 j = 1, 2*nlmax 
      a(i,j) = 0.0
23016 continue
23017 continue
23014 continue
23015 continue
      noff = 0
      i2flag = 0
      do23018 ilyr=1,nlyrs
      pertb(ilyr) = 1.01
      perta(ilyr) = 1.01
      pertr(ilyr) = 1.01
23018 continue
23019 continue
      do23020 iseis = 1, ns 
      loff = 0
      rp=p(iseis)
      dts=dt(iseis)
      dlys=dly(iseis)
      agauss=gauss(iseis)
      call partials( aa, rp, perta, pertb, pertr, nlyrs, nlmax, dts, ntm
     *ax, dlys, agauss, alpha, beta, rho, thk, pon )
      if(.not.(i2flag.eq.1))goto 23022
      return
23022 continue
      if(.not.(hpfilter))goto 23024
      do23026 j = 1, nlyrs-1 
      call hpbufilter(aa(1,j),nt(iseis),dt(iseis),fmin,npasses)
23026 continue
23027 continue
      call hpbufilter(aa(1,nlmax),nt(iseis),dt(iseis),fmin,npasses)
23024 continue
      do23028 i = 1, nt(iseis) 
      do23030 j = 1, nlyrs-1 
      a( i+noff, j+loff ) = aa( i, j )
23030 continue
23031 continue
      syn( i, iseis ) = aa( i, nlmax )
      b( i+noff ) = seis( i, iseis ) - syn( i, iseis )
      if(.not.( jumping ))goto 23032
      do23034 j = 1, nlyrs-1 
      b( i+noff ) = b( i+noff ) + aa(i,j) * beta(j)
23034 continue
23035 continue
23032 continue
23028 continue
23029 continue
      rms(iseis) = 0.0
      do23036 j = 1, nt(iseis) 
      rms(iseis) = rms(iseis) + ( seis(j,iseis) - syn(j,iseis) ) ** 2
23036 continue
23037 continue
      misfit(iseis) = rms(iseis) / sig_power(iseis)
      rms(iseis) = sqrt( rms(iseis) / nt(iseis) )
      noff = noff + nt(iseis)
23020 continue
23021 continue
      ruffa = 0.0
      ruffb = 0.0
      do23038 i = 1, nlyrs-2 
      ruffa = ruffa + ( alpha(i) - 2.*alpha(i+1) + alpha(i+2) ) ** 2
      ruffb = ruffb + ( beta(i) - 2.*beta(i+1) + beta(i+2) ) ** 2
23038 continue
23039 continue
      ruffa0 = sqrt( ruffa / (nlyrs-2) )
      ruffb0 = sqrt( ruffb / (nlyrs-2) )
      write(*,*) '******************************************************
     *****'
      write(*,*) 'Iteration: ', iter
      write(*,'(/)')
      write(*,*) 'initial fractional square misfit: ', ( misfit(i), i = 
     *1, ns )
      write(*,*) 'initial rms errors: ', ( rms(i), i = 1, ns )
      write(*,*) 'initial roughness alpha, beta:', ruffa0, ruffb0
      write(*,*) '******************************************************
     *****'
      if(.not.( jumping ))goto 23040
      do23042 i = 1, nlyrs 
      b( i+noff ) = 0.0
      b( i+noff + nlyrs-1 ) = 0.0
      if(.not.( i .lt. nlyrs-1 ))goto 23044
      a( i+noff, i ) = sigjmpb
      a( i+noff, i+1 ) = -2. * sigjmpb
      a( i+noff, i+2 ) = sigjmpb
      goto 23045
23044 continue
      if(.not.(i .eq. nlyrs-1))goto 23046
      a(i+noff, i) = sigjmpb
      a(i+noff, i+1) = -sigjmpb
      goto 23047
23046 continue
      a( i+noff, nlyrs ) = 1.0
      b( i+noff ) = beta(nlyrs)
23047 continue
23045 continue
23042 continue
23043 continue
      noff = noff + nlyrs
23040 continue
23048 if(.not.( .not. happy ))goto 23049
      call putsyn( ns, syn, ntmax, nsmax, dt, dura, dly, gauss, p, nt, p
     *orsv, iter, invnum )
      call wrtsoln( nlyrs, alpha, beta, rho, thk, iter ,invnum )
      npb = noff
      if(.not.( jumping ))goto 23050
      ip = nlyrs
      if(.not.(ns .eq. 1))goto 23052
      ip = nlyrs
23052 continue
      goto 23051
23050 continue
      ip = (nlyrs-1)
      if(.not.(ns .eq. 1))goto 23054
      ip = nlyrs-1
23054 continue
23051 continue
      call svdrs( a, ndat, npb, ip, b, ndat, 1, s)
      write(*,*) '******************************************************
     *****'
      write(*,*) 'Iteration: ', iter+1
      if(.not.( jumping ))goto 23056
      call jsoln( a, ndat, npb, ip, b, ndat, s, sol, tfraction)
      call putsvalues(s,ip,invnum)
      do23058 i = 1, nlyrs 
      beta(i) = sol(i)
      alpha(i) = beta(i)*vp_over_vs(i)
      rho(i) = 0.32 * alpha(i) + 0.77
23058 continue
23059 continue
23056 continue
      do23060 i = 1, nlyrs
      if(.not.(beta(i) .le. 0.0))goto 23062
      write(*,*)'Oops - negative velocities, quitting.'
      write(*,*)'Try increasing smoothing weight or'
      write(*,*)'decreasing initial-model perturbation size.'
      write(*,*)'Watch out for really slow near-surface'
      write(*,*)'layers and large initial model perturbations.'
      stop
23062 continue
23060 continue
23061 continue
      do23064 i = 1, ndat 
      do23066 j = 1, 2*nlmax 
      a(i,j) = 0.0
23066 continue
23067 continue
23064 continue
23065 continue
      noff = 0
      do23068 iseis = 1, ns 
      loff = 0
      rp=p(iseis)
      dts=dt(iseis)
      dlys=dly(iseis)
      agauss=gauss(iseis)
      call partials( aa, rp, perta, pertb, pertr, nlyrs, nlmax, dts, ntm
     *ax, dlys, agauss, alpha, beta, rho, thk, pon )
      if(.not.(hpfilter))goto 23070
      do23072 j = 1, nlyrs-1 
      call hpbufilter(aa(1,j),nt(iseis),dt(iseis),fmin,npasses)
23072 continue
23073 continue
      call hpbufilter(aa(1,nlmax),nt(iseis),dt(iseis),fmin,npasses)
23070 continue
      do23074 i = 1, nt(iseis) 
      do23076 j = 1, nlyrs-1 
      a( i+noff, j+loff ) = aa( i, j )
23076 continue
23077 continue
      syn( i, iseis ) = aa( i, nlmax )
      b( i+noff ) = seis( i, iseis ) - syn( i, iseis )
      if(.not.( jumping ))goto 23078
      do23080 j = 1, nlyrs-1 
      b( i+noff ) = b( i+noff ) + aa(i,j) * beta(j)
23080 continue
23081 continue
23078 continue
23074 continue
23075 continue
      rms(iseis) = 0.0
      do23082 j = 1, nt(iseis) 
      rms(iseis) = rms(iseis) + ( seis(j,iseis) - syn(j,iseis) ) ** 2
23082 continue
23083 continue
      misfit(iseis) = rms(iseis) / sig_power(iseis)
      rms(iseis) = sqrt( rms(iseis) / nt(iseis) )
      noff = noff + nt(iseis)
23068 continue
23069 continue
      if(.not.( jumping ))goto 23084
      do23086 i = 1, nlyrs 
      b( i+noff ) = 0.0
      b( i+noff + nlyrs-1 ) = 0.0
      if(.not.( i .lt. nlyrs-1 ))goto 23088
      a( i+noff, i ) = sigjmpb
      a( i+noff, i+1 ) = -2. * sigjmpb
      a( i+noff, i+2 ) = sigjmpb
      goto 23089
23088 continue
      if(.not.(i .eq. nlyrs-1))goto 23090
      a(i+noff, i) = sigjmpb
      a(i+noff, i+1) = -sigjmpb
      goto 23091
23090 continue
      a( i+noff, nlyrs ) = 1.0
      b( i+noff ) = beta(nlyrs)
23091 continue
23089 continue
23086 continue
23087 continue
      noff = noff + nlyrs
23084 continue
      ruffa = 0.0
      ruffb = 0.0
      do23092 i = 1, nlyrs-2 
      ruffa = ruffa + ( alpha(i) - 2.*alpha(i+1) + alpha(i+2) ) ** 2
      ruffb = ruffb + ( beta(i) - 2.* beta(i+1) + beta(i+2) ) ** 2
23092 continue
23093 continue
      ruffa = sqrt( ruffa / (nlyrs-2) )
      ruffb = sqrt( ruffb / (nlyrs-2) )
      write(*,*) 'fractional square misfit: ',( misfit(i), i = 1, ns )
      write(*,*) 'rms errors: ',( rms(i), i = 1, ns )
      write(*,*) 'roughness alpha, beta:',ruffa, ruffb
      if(.not.(ruffa0 .ne. 0 .and. ruffb0 .ne. 0))goto 23094
      write(*,'(a40,f8.2,1x,f8.2)') ' Percent Roughness Change (alpha,be
     *ta): ', 100*ruffa/ruffa0, 100*ruffb0/ruffb
      write(*,*) '******************************************************
     *****'
23094 continue
      iter = iter + 1
      happy = ( iter .ge. maxiter)
      goto 23048
23049 continue
      call wrtsoln( nlyrs, alpha, beta, rho, thk, iter, invnum )
      call putsyn( ns, syn, ntmax, nsmax, dt, dura, dly, gauss, p, nt, p
     *orsv, iter, invnum)
      write(*,*) 'Final Model '
      write(*,*) 'Layer   Vp/Vs       Poissons Ratio'
      do23096 ilyr=1,nlyrs
      vp_over_vs(ilyr) = alpha(ilyr)/beta(ilyr)
      pratio(ilyr) = vpovs_to_pr(vp_over_vs(ilyr))
      write(*,'(i5,2x,f10.7,8x,f5.3)') ilyr, vp_over_vs(ilyr),pratio(ily
     *r)
23096 continue
23097 continue
      return
      end
