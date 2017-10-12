      function ask(quest)
c
c   interactive i-o for real numbers
c
      character quest*(*)
      integer ounit
      common /innout/ inunit,ounit
      write(ounit,100) (quest(j:j),j=1,len(quest))
      read(inunit,*) anser
      ask=anser
      return
C*  100 format(80(a1,$))
  100 format(80(a1))
      end
