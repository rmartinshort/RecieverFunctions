      program pwavezrt
c  ****************************************************************
c  *    rotate zne components into zrt direction                  *
c  ****************************************************************
      complex data(4500,3)
      dimension caz(3)
      character eqfile*80,comp(3,2)*2
      character*32 evnam,stanam,outfil
      integer blank,ounit
      logical yes,yesno
      common /innout/ inunit,ounit
      data
     *     comp/'.z','.n','.e','.z','.r','.t'/
      inunit=5
      ounit=6
      pi=3.1415926
    2 call asktxt('specify event name: ',evnam)
      call asktxt('specify station name: ',stanam)
      e1=ask('pwave ending time after a (in Sec): ')
      iblank1=blank(evnam)
      iblank2=blank(stanam)
      iblank=iblank1+iblank2+25
      eqfile(1:22)='/var/data/field/pwave/'
      outfil(1:32)='                                '
      do 1 i=1,3
         eqfile(1:iblank)=eqfile(1:22)//evnam(1:iblank1)
     *   //'/'//stanam(1:iblank2)//comp(i,1)
         call zero(data(1,i),1,9000)
         call rsac1(eqfile,data(1,i),npts,beg,dt,9000,nerr)
         call getfhv('cmpaz',caz(i),nerr)
  1   continue
         call getfhv('a',a0,nerr)
         call getfhv('baz',baz,nerr)
       k11=(a0-10.-beg)/dt+1
       k12=(e1+a0-beg)/dt+1
       beg=-10.
       do 3 i=1,3
       call getdata(data(1,i),k11,k12,npts)
  3    continue
       npts=k12-k11+1
       call rotate(data,9000,3,baz,caz,npts)
       call minmax(data(1,1),npts,dmin,dmax,dmean)
       do 4 i=1,3
       do 5 j=1,npts
       data(j,i)=data(j,i)/dmax
  5    continue
  4    continue
       do 103 i=1,3
       outfil(1:iblank1+2)=evnam(1:iblank1)//comp(i,2)
       call wsac1(outfil,data(1,i),npts,beg,dt,nerr)
 103   continue
       yes=yesno('Try another (y or n)? ')
       if(yes) go to 2
       end
