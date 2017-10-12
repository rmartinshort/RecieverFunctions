      subroutine crossfun(x,y,cross,n1,n2)
      dimension x(1),y(1),cross(n1),x1(8000)
      do 20 i=1,n1+n2
      x1(i)=0.0
 20   continue
      do 30 i=1,n2
      x1(i)=y(i)
 30   continue
      do 10 k=1,n1
      cross(k)=0.0
      do 10 i=1,n1
      cross(k)=cross(k)+x(i)*x1(i+k-1)
 10   continue
      return
      end
