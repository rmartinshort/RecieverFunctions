        dimension data(8000),a(8000),v(8000),c(8000),r(8000)
        character eqfile*32
        integer ounit
        common /innout/inunit,ounit
        inunit=5
        ounit=6
        call asktxt('specify input file : ',eqfile)
        call rsac1(eqfile,data,npts,beg,dt,8000,nerr)
        call crossfun1(data,data,r,npts,npts)
        do 1 i=1,npts
        r(i)=r(i)/r(1)
 1      continue
        c(1)=-1.
        a(1)=1.
        v(1)=1.
        do 220 j=2,npts
        a(j)=0.
        e=0.
        do 210 i=2,j
 210    e=e+r(i)*a(j-i+1)
        c(j)=e/v(j-1)
        v(j)=v(j-1)-e*c(j)
        jh=(j+1)/2
        do 220 i=1,jh
        bot=a(j-i+1)-c(j)*a(i)
        a(i)=a(i)-c(j)*a(j-i+1)
 220    a(j-i+1)=bot
        do 2 i=2,npts
        a(i-1)=a(i)
 2      continue
        call wsac1('xxx',a,npts,beg,dt,nerr)
        end
