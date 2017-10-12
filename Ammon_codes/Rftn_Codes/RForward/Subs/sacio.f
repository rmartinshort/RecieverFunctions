      subroutine sacio(file,x,np,dt,inout)
      include './mach'
      include './hdr'
      integer year,jday,hour,min,isec,msec,ounit
      character*8 sta,compnm,evnm
      character file*64
      dimension x(1)
      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *              compnm,caz,cinc,evnm,bz,del,rayp,dep,decon,agauss,
     *              c,tq,rinstr,dlen,btimey,ty0,ty1,ty2
      common /innout/ inunit,ounit
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c
c **************************************************************************
c
c   parameters are:
c
c      dmin,dmax,dmean = min,max, and mean of data read or written
c      year,jday,hour,min,isec,msec = gmt reference time (all integers)
c      sta = station name (a8)
c      compnm = component name (a8)
c      caz = orientation of component (wrt north)
c      cinc = inclination of component (wrt vertical)
c      evnm = name of event (a8)
c      bz  = back azimuth of from station to event
c      del = distance from station to event in degrees
c      rayp  = ray parameter of arriving phase in j-b earth
c      dep = depth of event in kilometers
c      decon = if = 1. indicates data has been source equalized
c      agauss = width of gaussian used in source equalization if decon = 1.
c      c     = trough filler used in source equalization if decon = 1.
c      tq = t/q value used in synthetic (if used)
c      rinstr = 1. if response of 15-100 system has been put into synthetic
c      dlen  = length of data in secs.
c      btimey = time 0f 1st data point wrt gmt reference time (in secs)
c      ty0,ty1,ty2 = user defined times wrt gmt reference (in secs)
c
c ****************************************************************************
c
c    call to sacio is:
c                      call sacio(file,x,np,dt,inout)
c
c    where file = file to be read or written
c          x    = data array to be used
c          np   = number of points in x
c          dt   = sampling rate for x
c          inout = +1 for reading a sac file
c                = -1 for writing a sac file
c
c ******************************************************************************
      if(inout.lt.0) goto 1
c
c  read a sac file
c
c     call inicm
c     call inihdr
      call rsac1(file,x,np,b,dt,100000,nerr)
      if(nerr.eq.0) go to 2
      write(ounit,100) nerr
  100 format(' nerr = ',i2,' in sacio read')
      return
2     continue

      call getkhv('kstnm',sta,nerr)
      call getfhv('depmin',dmin,nerr)
      call getfhv('depmax',dmax,nerr)
      call getfhv('depmen',dmean,nerr)

      call getnhv('nzyear',year,nerr)
      call getnhv('nzhour',hour,nerr)
      call getnhv('nzjday',jday,nerr)
      call getnhv('nzmin',  min,nerr)
      call getnhv('nzsec', isec,nerr)
      call getnhv('nzmsec',msec,nerr)

      btimey=b
      call getfhv('e',e,nerr)
      dlen=e-b

      call getkhv('kcmpnm',compnm,nerr)
      call getfhv('cmpaz',caz,nerr)
      call getfhv('cmpinc',cinc,nerr)
      call getkhv('kcmpnm',compnm,nerr)
      call getfhv('baz',bz,nerr)
      call getfhv('gcarc',del,nerr)
      call getfhv('user0',rayp,nerr)
      call getfhv('user1',dep,nerr)
      call getfhv('user2',agauss,nerr)
      call getfhv('user3',c,nerr)
      if(c.gt.0.) decon=1.0
      call getfhv('user4',tq,nerr)
      call getfhv('user5',rinstr,nerr)
      ty0=t0
      ty1=t1
      ty2=t2
      return
c
c write a sac file
c
    1 continue
      call newhdr
      call setnhv('npts',np,nerr)
      call setkhv('kstnm',sta,nerr)
      call setfhv('delta',dt,nerr)
      call setfhv('depmin',dmin,nerr)
      call setfhv('depmax',dmax,nerr)
      call setfhv('depmen',dmean,nerr)
      if(year.lt.1960) go to 4
      if(decon.lt..0001) go to 5
      call setfhv('b',btimey,nerr)
      go to 6
    5 continue
      call setnhv('nzyear',year,nerr)
      call setnhv('nzhour',hour,nerr)
      call setnhv('nzjday',jday,nerr)
      call setnhv('nzmin',min,nerr)
      call setnhv('nzsec',isec,nerr)
      call setnhv('nzmsec',msec,nerr)
    4 continue
      call setfhv('b',btimey,nerr)

    6 continue
      if(abs(ty0).gt..00001) t0=ty0
      if(abs(ty1).gt..00001) t1=ty1
      if(abs(ty2).gt..00001) t2=ty2
      call setkhv('kcmpnm',compnm,nerr)
      call setfhv('cmpaz',caz,nerr)
      call setfhv('cmpinc',cinc,nerr)
      call setkhv('kevnm',evnm,nerr)
      call setfhv('baz',bz,nerr)
      if(del.gt..0001) call setfhv('gcarc',del,nerr)
      call setfhv('user0',rayp,nerr)
      if(dep.gt..0001) call setfhv('user1',dep,nerr)
      if(decon.lt..0001) go to 3
      call setfhv('user2',agauss,nerr)
      call setfhv('user3',c,nerr)
    3 if(tq.gt..0001) call setfhv('user4',tq,nerr)
      if(rinstr.gt..0001) call setfhv('user5',rinstr,nerr)
      call wsac0(file,xxx,x,nerr)
      if(nerr.eq.0) return
      write(ounit,101) nerr
  101 format(' nerr = ',i2,' in sacio write')
      return
      end
