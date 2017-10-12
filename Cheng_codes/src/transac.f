       subroutine transac(infil,outfil,b0,b1,x)
       dimension x(1)
       character infil*32,outfil*32
       include 'hdr.inc'
       inunit1=8
       open(inunit1,file=infil,access='direct',recl=158*4)
       read(inunit1,rec=1)(rheader(i),i=1,70),(iheader(i),i=1,35),
     *  (lheader(i),i=1,5),kstnm,kevnm,(cheader(i),i=1,21)
       close(inunit1)
c       n1=(rheader(9)-rheader(6)-b0)/rheader(1)
       b1=rheader(7)
       n1=1
c       n2=(rheader(9)-rheader(6)+b1)/rheader(1)+1
       n2=(rheader(7)-rheader(6))/rheader(1)+1
       open(inunit1,file=infil,access='direct',recl=(158+n2)*4)
       read(inunit1,rec=1)(idummy,i=1,158+n1-1),(x(i-n1+1),i=n1,n2)
       close(inunit1)
       b0=rheader(6)-rheader(9)
       rheader(6)=-1*b0
       rheader(7)=b1
       rheader(9)=0.
       iheader(10)=n2-n1+1
       do 1 i=1,6
       iheader(i)=-12345
   1   continue
       call minmax(x,iheader(10),rheader(2),rheader(3),rheader(57))
       open(inunit1,file=outfil,access='direct',
     *   recl=(158+iheader(10))*4)
       write(inunit1,rec=1)(rheader(i),i=1,70),(iheader(i),i=1,35),
     *  (lheader(i),i=1,5),kstnm,kevnm,(cheader(i),i=1,21),
     *  (x(i),i=1,iheader(10))
       close(inunit1)
       end
