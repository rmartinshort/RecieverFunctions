      subroutine sacio(file,x,np,dt,inout)
      integer inunit,ounit
      integer year,jday,hour,min,isec,msec
      character*8 sta,compnm,evnm
      common /tjocm/dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *             compnm,caz,cinc,evnm,bz,del,rayp,dep,decon,agauss,
     *              c,tq,rinstr,dlen,btimey,ty0,ty1,ty2
      character file*32
      dimension x(1)
      common /innout/ inunit,ounit
      include 'hdr.inc'
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
      call rsac1(file,x,np,btimey,dt,10000,nerr)
      if(nerr.eq.0) go to 2
      write(ounit,100) nerr
  100 format('nerr= ',i2,'in sacio read')
    2 continue
      call rhdr(file,np,dt)
      return
c
c write a sac file
C
c    1 call wsac0(file,xxx,x,nerr)
    1 call wsac1(file,x,np,btimey,dt,nerr)
      if(nerr.eq.0) go to 3
      write(ounit,200) nerr
  200 format('nerr=',i2,'in sacio write')
    3 continue
      call whdr(file,np,dt)
      return
      end
