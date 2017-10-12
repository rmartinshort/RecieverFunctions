      integer year,jday,hour,min,isec,msec
      dimension x(6000)
      character*8 sta,compnm,evnm
      character file1*32
      common /tjocm/ dmin,dmax,dmean,year,jday,hour,min,isec,msec,sta,
     *              compnm,caz,cinc,evnm,bz,del,rayp,dep,decon,agauss,
     *              c,tq,rinstr,dlen,btimey,ty0,ty1,ty2,np,dt,fia
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
        file1='amdo.e'
        call rhdr(file1,x,80.,150.)
c        call rdsac(file1, x,np)
c        write(*,*) (rheader(i),i=1,70)
c        write(*,*) (iheader(i),i=1,35)
c        write(*,*) (lheader(i),i=1,5)
c        write(*,*) kstnm,kevnm
c        write(*,*) (cheader(i),i=1,21)
        file1='wndo1.z'
c        write(*,*) 'np=',np,'dt=',dt
        call whdr
        call wdsac(file1,x,np)
        kevnm='b'
        kstnm='c'
        cheader(3)='cccccccc'
        e=btimey+dlen
        call rhdr(file1,x,btimey,e)
        write(*,*) (rheader(i),i=1,70)
        write(*,*) (iheader(i),i=1,35)
        write(*,*) (lheader(i),i=1,5)
        write(*,*)  kstnm,kevnm
        write(*,*) (cheader(i),i=1,21)
cc        write(*,*) (x(i),i=1,np)
        end
