c      subroutine csoln( a, NDAT, npb, ip, b, NDAT, s, sol )
      subroutine csoln( a, NDAT, npb, ip, b,  s, sol )
      integer NDAT
      integer npb,ip
      real a(NDAT,ip), b(NDAT), s(ip), sol(ip)

c-- initialize soln
      do 10 i = 1, ip 
        sol(i) = 0.0
  10  continue
      do 20 j = 1, ip 
C-- if the singular value is zero, drop the solution
c   #-- (use the ordered eigen values to break??
        if ( s(j) .ne. 0.0 )  then
           p = b(j) / s(j)
        else 
           p = 0.0
        endif
        do 30 i = 1, ip 
           sol(i) = sol(i) + p * a(i,j)
  30    continue
  20  continue
      return
      end
