       dimension x(10000)
       call rsac1('021jz_00106.bhe',x(1),npts,beg,dt,9000,nerr)
       write(*,*) nerr, npts, beg, dt
       end

