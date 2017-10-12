      parameter (n=2)
      real twobytwo(2,2) /4*-1/
      call mkidentity (twobytwo,n)
      write(*,*) 'Ihere'
      print *, determinant(twobytwo)
      end
