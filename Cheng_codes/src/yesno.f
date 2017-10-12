      logical function yesno(quest)
c
c   interactive i-o for logical variables
c    yesno must be declared logical in calling program
c
      character quest*(*),answer*1
      logical lanswr
      integer ounit
      common /innout/ inunit,ounit
      write(ounit,100) (quest(j:j),j=1,len(quest))
      read(inunit,200) answer
      lanswr=.false.
      if(answer.eq.'y') lanswr=.true.
      yesno=lanswr
C*  100 format(80(a1,$))
  100 format(80(a1))
  200 format(a1)
      return
      end
