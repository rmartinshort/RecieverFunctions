      subroutine crossfun(x,y,cross,n1,n2)
      dimension x(1),y(1),x1(8000),cross(n1+n2)
      do 1 i=1,n2+2*n1
      x1(i)=0.0
      cross(i)=0.0
 1    continue
      do 2 i=n1+1,n1+n2
      x1(i)=y(i-n1)
 2    continue
      do 10 k=1,n1+n2
      cross(k)=0.0
      do 10 i=1,n1
      cross(k)=cross(k)+x(i)*x1(i+k)
 10   continue
      return
      end
