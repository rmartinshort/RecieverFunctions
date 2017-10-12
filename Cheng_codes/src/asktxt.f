      subroutine asktxt(quest,answer)
c
c   interactive i-o for character strings
c      string returned may be a maximun of 32 characters
c
      character answer*32,quest*(*)
      integer ounit
      common /innout/ inunit,ounit
      write(ounit,100) (quest(j:j),j=1,len(quest))
      read(inunit,200) answer
C*  100 format(80(a1,$))
  100 format(80(a1))
  200 format(a32)
      return
      end
