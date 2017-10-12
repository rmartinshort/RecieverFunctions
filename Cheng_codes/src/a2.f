      subroutine mkidentity (array,m)
      real array(m,m)
      do 10 i=1,m
      do 20 j=1,m
      if (i.eq.j) then
         array(i,j)=1.
      else
         array(i,j)=0.
      endif
 20   continue
 10   continue
      return
      end
