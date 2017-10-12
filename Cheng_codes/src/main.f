         dimension a(5,5),v(5,5),b(5),z(5),d(5)
         logical ev
         write(*,*) 'input ev:'
         read(*,*) ev
         write(*,*) 'ev=',ev
         n=5
         a(1,1)=10
         a(1,2)=1
         a(1,3)=2
         a(1,4)=3
         a(1,5)=4
         a(2,1)=1
         a(2,2)=9
         a(2,3)=-1
         a(2,4)=2
         a(2,5)=-3
         a(3,1)=2
         a(3,2)=-1
         a(3,3)=7 
         a(3,4)=3
         a(3,5)=-5
         a(4,1)=3
         a(4,2)=2
         a(4,3)=3
         a(4,4)=12
         a(4,5)=-1
         a(5,1)=4
         a(5,2)=-3
         a(5,3)=-5
         a(5,4)=-1
         a(5,5)=15
         do 12 i=1,n
         do 11 j=1,n
         z(j)=a(i,j)
 11      continue
         write(*,*) z
 12      continue
         call jacobi(n,ev,a,d,v,b,z,kk)
         write(*,*) 'kk=',kk
c         do 20 i=1,5
c         write(*,*) 'i=',i
c         write(*,*) (a(i,j),j=1,5)
c  20     continue
         write(*,*) 'd=',d
         do 30 j=1,5
         write(*,*) 'j=',j
         write(*,*) (v(i,j),i=1,5)
  30     continue
         end
         
