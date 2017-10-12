      program predecon
c
c  ********************************************************************
c
c    fortran 77 program to perpare for  source equalization deconvolution
c    on a three component seismogram,fistly,obtain the pwaveform,secondly
c    ,rotate the horizontal components into backazimuth 
c
c  *************************************************************************
c
      dimension caz(3)
      character infil*32,outfil*32,rotfil*32,comp(3,2)*6
      dimension data(3000,3)
      logical yes,yesno,rel
      integer blank,ounit
c **********************************************************************
      include 'hdr.inc'
c
c **************************************************************************
c
c   parameter definitions may be found in sac comments
c
c      common /win/ data
      common /innout/ inunit,ounit
      data
     *     comp/'.z    ','.n    ','.e    ','.z    ','.r    ','.t    '/
      inunit=5
      ounit=6
      outfil='                                '
      rotfil='                                '
    2 call asktxt('Specify transform file: ',infil)
      rel=yesno('Geographic(zne) data (y or n)? ')
      iblank=blank(infil)
      call asktxt('specify outfil: ',outfil)
      iblank1=blank(outfil)
      isyntp=2
      if(rel) isyntp=1
      write(*,*) 'input end time :'
      read(*,*) b1
      do 1 i=1,3
         call zero(data(1,i),1,3000)
         infil(1:iblank + 6)=infil(1:iblank)//comp(i,isyntp)
         outfil(1:iblank1+7)=outfil(1:iblank1)//'p'//comp(i,isyntp)
         call transac(infil,outfil,b0,b1,data(1,i))
         caz(i)=rheader(58)
    1 continue
      npts=iheader(10)
      if(rel) call rotate(data,3000,3,rheader(53),caz,npts)
      do 5 i=1,3
      rotfil(1:iblank1+7)=outfil(1:iblank1)//'r'//comp(i,2)
      if(i.eq.2) rheader(58)=rheader(53)+180
      if(i.eq.3) rheader(58)=rheader(53)+270
      if(rheader(58).gt.360.) rheader(58)=rheader(58)-360.
      rheader(59)=90.
      call minmax(data(1,i),npts,rheader(2),rheader(3),rheader(57))
      call wdsac(rotfil,data(1,i))
    5 continue
c
c    *******************************************************************
c
      yes=yesno('Try another (y or n)? ')
      if(yes) go to 2
      stop
      end
