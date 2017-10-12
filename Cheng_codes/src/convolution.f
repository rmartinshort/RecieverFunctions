       subroutine  convolution(x,y,z,n1,n2)
c  ****************************************************************
c  *     This is the convolution program between x(n1) and y(n2)  *
c  *     the output is z(n1+n2-1)                                 *
c  *                 n1    n2                                     *
c  *     z(i+j-1) = sum ( sum (x(i)*y(j)) )                       *
c  *                 i=1  j=1                                     *
c  ****************************************************************
       dimension x(1),y(1),z(n1+n2-1)
c       do 1 i=1,n2
c       z(i)=0.
c       do 1 j=1,i
c       z(i)=z(i)+x(j)*y(i-j+1)
c 1     continue
c       do 2 i=2,n1
c       k=n2+i-1
c       z(k)=0.
c       do 2 j=i,n1
c       z(k)=z(k)+x(j)*y(n2+i-j)
c  2    continue
c       k0=n1+n2-1
       do 10 i=1,n1+n2-1
       z(i)=0.
  10   continue
       do 20 i=1,n1
       do 20 j=1,n2
       z(i+j-1)=z(i+j-1)+x(i)*y(j)
  20   continue
       return
       end
