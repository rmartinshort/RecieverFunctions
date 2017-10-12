
c   imsl routine name   - lsvdf
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - june 1, 1980
c
c   purpose             - singular value decomposition of a real
c                           matrix
c
c   usage               - call lsvdf (a,ia,m,n,b,ib,nb,s,wk,ier)
c
c   arguments    a      - real m by n matrix. (input/output)
c                         on input, a contains the matrix to be
c                           decomposed.
c                         on output, a contains the n by n matrix v
c                           in its first n rows. see remarks. either
c                           m.ge.n or m.lt.n is permitted.
c                ia     - row dimension of matrix a exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                           a is used by lsvdf as work storage for an
c                           n by n matrix. therefore, ia must be
c                           greater than or equal to max(m,n).
c                m      - number of rows in a. (input)
c                n      - number of columns in a. (input)
c                b      - m by nb matrix. (input/output)
c                           b is not used if nb.le.0. otherwise, b is
c                           replaced by the matrix product u**(t) * b.
c                           see remarks.
c                ib     - row dimension of matrix b exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                nb     - number of columns in b. (input)
c                           if nb.le.0, b is not used.
c                s      - vector of length n. (output)
c                         on output, s contains the ordered singular
c                           values of a.  s(1) .ge. s(2),...,
c                           .ge. s(n) .ge. 0.
c                wk     - work vector of length 2n.
c                ier    - error parameter. (output)
c                         warning error
c                           ier=33 indicates that matrix a is not
c                             full rank or very ill-conditioned. small
c                             singular values may not be very accurate.
c                           ier=34 indicates that either n.le.0 or
c                             m.le.0.
c                         terminal error
c                           ier=129 indicates that convergence was
c                             not obtained by lsvdb and computation
c                             was discontinued.
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - lsvdb,lsvg1,lsvg2,vhs12,uerset,uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks  1.  lsvdf computes the singular value decomposition of
c                a real m by n matrix
c                     a = u * q * v**(t)  where
c                u is an m by m orthogonal matrix,
c                v is an n by n orthogonal matrix, and
c                q is an m by n matrix with all elements zero except
c                     q(i,i) = s(i) i=1,...,min(m,n).
c                v is returned in the first n rows of a.
c                u is obtained by setting b to the m by m identity
c                matrix, on input, and setting nb=m. on output, b is
c                replaced by u**(t).
c            2.  the notation u**(t) and v**(t) represents u
c                transpose and v transpose, respectively. q**(+)
c                denotes the generalized inverse of q.
c            3.  lsvdf is useful in analyzing and solving the least
c                squares problem a*x.appr.b (where .appr. means
c                approximately equals). in this case b is a vector of
c                length m and lsvdf is called with ib=m, nb=1. the
c                solution is x=v*q**(+)*u**(t)*b. u**(t)*b replaces
c                b on output. the solution x is obtained as follows...
c                (the user may wish to set small singular values such as
c                s(i) to zero just prior to computing x. see reference
c                for details.)
c
c                c                  compute q**(+) * u**(t) * b
c                      l=min0(m,n)
c                      do 10 i=1,l
c                         t=0.0
c                         if (s(i).ne.0.0) t=b(i)/s(i)
c                         b(i)=t
c                   10 continue
c                c                  compute v * q**(+) * u**(t) * b
c                   15 do 25 i=1,n
c                         x(i)=0.0
c                         do 20 j=1,l
c                   20    x(i)=x(i)+a(i,j)*b(j)
c                   25 continue
c                if b is set to the j-th column of the m by m identity
c                matrix on input, x is the j=th column of a**(+) (the
c                generalized inverse of a).
c            4.  the user should be aware of several practical aspects
c                of the singular value analysis. one of these is the
c                effect of the uncertainty of the data and the need for
c                scaling. see the lawson-hanson reference pages 180-198
c                for details.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine lsvdf  (a,ia,m,n,b,ib,nb,s,wk,ier)
c
c                                  specifications for arguments
      integer            ia,m,n,ib,nb,ier
      real               a(ia,n),b(ib,1),s(n),wk(n,2)
c                                  specifications for local variables
      integer            i,j,jp1,k,l,mm,nn,nnp1,ns,nsp1
      real               zero,one,t
      data               zero/0.0/,one/1.0/
c                                  first executable statement
      ier = 0
c                                  begin special for zero rows and
c                                    cols. pack the nonzero cols to the
c                                    left
      nn = n
      ier = 34
      if (nn.le.0.or.m.le.0) go to 9000
      ier = 0
      j = nn
    5 continue
      do 10 i=1,m
         if (a(i,j).ne.zero) go to 25
   10 continue
c                                  col j is zero. exchange it with col
c                                    n
      if (j.eq.nn) go to 20
      do 15 i=1,m
   15 a(i,j) = a(i,nn)
   20 continue
      a(1,nn) = j
      nn = nn-1
   25 continue
      j = j-1
      if (j.ge.1) go to 5
c                                  if n=0 then a is entirely zero and
c                                    svd computation can be skipped
      ns = 0
      if (nn.eq.0) go to 120
c                                  pack nonzero rows to the top quit
c                                    packing if find n nonzero rows
      i = 1
      mm = m
   30 if (i.gt.n.or.i.ge.mm) go to 75
      if (a(i,i).ne.zero) go to 40
      do 35 j=1,nn
         if (a(i,j).ne.zero) go to 40
   35 continue
      go to 45
   40 i = i+1
      go to 30
c                                  row i is zero exchange rows i and m
   45 if (nb.le.0) go to 55
      do 50 j=1,nb
         t = b(i,j)
         b(i,j) = b(mm,j)
         b(mm,j) = t
   50 continue
   55 do 60 j=1,nn
   60 a(i,j) = a(mm,j)
      if (mm.gt.nn) go to 70
      do 65 j=1,nn
   65 a(mm,j) = zero
   70 continue
c                                  exchange is finished
      mm = mm-1
      go to 30
c
   75 continue
c                                  end special for zero rows and
c                                    columns
c                                  begin svd algorithm..
c                                  (1) reduce the matrix to upper
c                                    bidiagonal form with householder
c                                    transformations.
c                                    h(n)...h(1)aq(1)...q(n-2) =
c                                    (d**t,0)**t where d is upper
c                                    bidiagonal.
c                                  (2) apply h(n)...h(1) to b. here
c                                    h(n)...h(1)*b replaces b in
c                                    storage.
c                                  (3) the matrix product w=
c                                    q(1)...q(n-2) overwrites the first
c                                    n rows of a in storage.
c                                  (4) an svd for d is computed. here k
c                                    rotations ri and pi are computed
c                                    so that rk...r1*d*p1**(t)...pk**(t)
c                                    = diag(s1,...,sm) to working
c                                    accuracy. the si are nonnegative
c                                    and nonincreasing. here rk...r1*b
c                                    overwrites b in storage while
c                                    a*p1**(t)...pk**(t) overwrites a
c                                    in storage.
c                                  (5) it follows that,with the proper
c                                    definitions, u**(t)*b overwrites
c                                    b, while v overwrites the first n
c                                    row and columns of a.
      l = min0(mm,nn)
c                                  the following loop reduces a to
c                                    upper bidiagonal and also applies
c                                    the premultiplying transformations
c                                    to b.
      do 85 j=1,l
         if (j.ge.mm) go to 80
         jp1 = min0(j+1,nn)
         call vhs12 (1,j,j+1,mm,a(1,j),1,t,a(1,jp1),1,ia,nn-j)
         call vhs12 (2,j,j+1,mm,a(1,j),1,t,b,1,ib,nb)
   80    if (j.ge.nn-1) go to 85
         call vhs12 (1,j+1,j+2,nn,a(j,1),ia,wk(j,2),a(j+1,1),ia,1,mm-j)
   85 continue
c                                  copy the bidiagonal matrix into the
c                                    array s for lsvdb
      if (l.eq.1) go to 95
      do 90 j=2,l
         s(j) = a(j,j)
         wk(j,1) = a(j-1,j)
   90 continue
   95 s(1) = a(1,1)
c
      ns = nn
      if (mm.ge.nn) go to 100
      ns = mm+1
      s(ns) = zero
      wk(ns,1) = a(mm,mm+1)
  100 continue
c                                  construct the explicit n by n
c                                    product matrix, w=q1*q2*...*ql*i
c                                    in the array a
      do 115 k=1,nn
         i = nn+1-k
         if (i.gt.min0(mm,nn-2)) go to 105
         call vhs12 (2,i+1,i+2,nn,a(i,1),ia,wk(i,2),a(1,i+1),1,ia,nn-i)
  105    do 110 j=1,nn
  110    a(i,j) = zero
         a(i,i) = one
  115 continue
c                                  compute the svd of the bidiagonal
c                                    matrix
c
      level=1
      call uerset(level,levold)
      call lsvdb (s(1),wk(1,1),ns,a,ia,nn,b,ib,nb,ier)
c                                  test for ier=33
c
      if (ier.gt.128) go to 9000
      call uerset(levold,levold)
      if (ier.ne.33) go to 120
      t=0.0
      nm=min0(m,n)
      if (s(1).ne.zero) t=s(nm)/s(1)
      f=100.0+t
      if (f.eq.100.0) go to 120
      ier=0
  120 continue
      if (ns.ge.nn) go to 130
      nsp1 = ns+1
      do 125 j=nsp1,nn
  125 s(j) = zero
  130 continue
      if (nn.eq.n) go to 155
      nnp1 = nn+1
c                                  move record of permutations and
c                                    store zeros
      do 140 j=nnp1,n
         s(j) = a(1,j)
         if (nn.lt.1) go to 140
         do 135 i=1,nn
  135    a(i,j) = zero
  140 continue
c                                  permute rows and set zero singular
c                                    values
      do 150 k=nnp1,n
         i = s(k)
         s(k) = zero
         do 145 j=1,n
            a(k,j) = a(i,j)
  145    a(i,j) = zero
         a(i,k) = one
  150 continue
c                                  end special for zero rows and
c                                    columns
  155 if (ier.eq.0) go to 9005
 9000 continue
      call uertst (ier,6hlsvdf )
 9005 return
      end
c   imsl routine name   - lsvdb
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - singular value decomposition of a bidiagonal
c                           matrix.
c
c   usage               - call lsvdb (d,e,n,v,iv,nrv,c,ic,ncc,ier)
c
c   arguments    d      - vector of length n. (input/output)
c                         on input, d contains the diagonal elements
c                           of the bidiagonal matrix b. d(i)=b(i,i),
c                           i=1,...,n.
c                         on output, d contains the n (nonnegative)
c                           singular values of b in nonincreasing
c                           order.
c                e      - vector of length n. (input/output)
c                         on input, e contains the superdiagonal
c                           elements of b. e(1) is arbitrary,
c                           e(i)=b(i-1,i), i=2,...,n.
c                         on output, the contents of e are modified
c                           by the subroutine.
c                n      - order of the matrix b. (input)
c                v      - nrv by n matrix. (input/output)
c                           if nrv.le.0, v is not used. otherwise,
c                           v is replaced by the nrv by n product
c                           matrix v*vb. see remarks.
c                iv     - row dimension of matrix v exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                nrv    - number of rows of v. (input)
c                c      - n by ncc matrix. (input/output)
c                           if ncc.le.0 c is not used. otherwise, c
c                           is replaced by the n by ncc product
c                           matrix ub**(t) * c. see remarks.
c                ic     - row dimension of matrix c exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c                ncc    - number of columns in c. (input)
c                ier    - error parameter. (input)
c                         warning error
c                           ier=33 indicates that matrix b is not full
c                             rank or very ill-conditioned. small
c                             singular values may not be very accurate.
c                         terminal error
c                           ier=129 indicates that convergence was
c                             not attained after 10*n qr sweeps.
c                             (convergence usually occurs in about
c                             2*n sweeps).
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - lsvg1,lsvg2,vhs12,uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks      lsvdb computes the singular value decomposition of
c                an n by n bidiagonal matrix
c                     b = ub * s * vb**(t)    where
c                ub and vb are n by n orthogonal matrices and
c                s is diagonal.
c                if arguments v and c are n by n identity matrices,
c                on exit they are replaced by vb and ub**t,
c                respectively.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine lsvdb (d,e,n,v,iv,nrv,c,ic,ncc,ier)
c
c                                  specifications for arguments
      integer            n,iv,nrv,ic,ncc,ier
      real               d(n),e(n),v(iv,1),c(ic,1)
c                                  specifications for local variables
      integer            i,ii,j,k,kk,l,ll,lp1,nqrs,n10
      logical            wntv,havers,fail
      real               dnorm,zero,one,two,cs,f,g,h,sn,t,x,y,z
      data               zero/0.0/,one/1.0/,two/2.0/
c                                  first executable statement
      ier = 0
      if (n.le.0) go to 9005
      n10 = 10*n
      wntv = nrv.gt.0
      havers = ncc.gt.0
      fail = .false.
      nqrs = 0
      e(1) = zero
      dnorm = zero
      do 5 j=1,n
    5 dnorm = amax1(abs(d(j))+abs(e(j)),dnorm)
      do 100 kk=1,n
         k = n+1-kk
c                                  test for splitting or rank
c                                    deficiencies first make test for
c                                    last diagonal term, d(k), being
c                                    small.
   10    if (k.eq.1) go to 25
         t = dnorm+d(k)
         if (t.ne.dnorm) go to 25
c
c                                  since d(k) is small we will make a
c                                    special pass to transform e(k) to
c                                    zero.
         cs = zero
         sn = -one
         do 20 ii=2,k
            i = k+1-ii
            f = -sn*e(i+1)
            e(i+1) = cs*e(i+1)
            t = d(i)
            call lsvg1 (t,f,cs,sn,d(i))
c                                  transformation constructed to zero
c                                    position (i,k).
            if (.not.wntv) go to 20
            do 15 j=1,nrv
   15       call lsvg2 (cs,sn,v(j,i),v(j,k))
c
c                                  accumulate rt. transformations in v.
   20    continue
c                                  the matrix is now bidiagonal, and of
c                                    lower order since e(k) .eq. zero
   25    do 30 ll=1,k
            l = k+1-ll
            t = dnorm+e(l)
            if (t.eq.dnorm) go to 50
            t = dnorm+d(l-1)
            if (t.eq.dnorm) go to 35
   30    continue
c                                  this loop cant complete since e(1) =
c                                    zero.
         go to 50
c                                  cancellation of e(l), l.gt.1.
   35    cs = zero
         sn = -one
         do 45 i=l,k
            f = -sn*e(i)
            e(i) = cs*e(i)
            t = dnorm+f
            if (t.eq.dnorm) go to 50
            t = d(i)
            call lsvg1 (t,f,cs,sn,d(i))
            if (.not.havers) go to 45
            do 40 j=1,ncc
   40       call lsvg2 (cs,sn,c(i,j),c(l-1,j))
   45    continue
c                                  test for convergence
   50    z = d(k)
         if (l.eq.k) go to 85
c                                  shift from bottom 2 by 2 minor of
c                                    b**(t)*b.
         x = d(l)
         y = d(k-1)
         g = e(k-1)
         h = e(k)
         f = ((y-z)*(y+z)+(g-h)*(g+h))/(two*h*y)
         g = sqrt(one+f**2)
         if (f.lt.zero) go to 55
         t = f+g
         go to 60
   55    t = f-g
   60    f = ((x-z)*(x+z)+h*(y/t-h))/x
c                                  next qr sweep
         cs = one
         sn = one
         lp1 = l+1
         do 80 i=lp1,k
            g = e(i)
            y = d(i)
            h = sn*g
            g = cs*g
            call lsvg1 (f,h,cs,sn,e(i-1))
            f = x*cs+g*sn
            g = -x*sn+g*cs
            h = y*sn
            y = y*cs
            if (.not.wntv) go to 70
c                                  accumulate rotations (from the
c                                    right) in v
            do 65 j=1,nrv
   65       call lsvg2 (cs,sn,v(j,i-1),v(j,i))
   70       call lsvg1 (f,h,cs,sn,d(i-1))
            f = cs*g+sn*y
            x = -sn*g+cs*y
            if (.not.havers) go to 80
            do 75 j=1,ncc
   75       call lsvg2 (cs,sn,c(i-1,j),c(i,j))
c
c                                  apply rotations from the left to
c                                    right hand sides in c
   80    continue
         e(l) = zero
         e(k) = f
         d(k) = x
         nqrs = nqrs+1
         if (nqrs.le.n10) go to 10
c                                  return to test for splitting.
         fail = .true.
c                                  cutoff for convergence failure. nqrs
c                                    will be 2*n usually.
   85    if (z.ge.zero) go to 95
         d(k) = -z
         if (.not.wntv) go to 95
         do 90 j=1,nrv
   90    v(j,k) = -v(j,k)
   95    continue
c                                  convergence. d(k) is made
c                                    nonnegative
  100 continue
      if (n.eq.1) go to 140
      do 105 i=2,n
         if (d(i).gt.d(i-1)) go to 110
  105 continue
      go to 140
c                                  every singular value is in order
  110 do 135 i=2,n
         t = d(i-1)
         k = i-1
         do 115 j=i,n
            if (t.ge.d(j)) go to 115
            t = d(j)
            k = j
  115    continue
         if (k.eq.i-1) go to 135
         d(k) = d(i-1)
         d(i-1) = t
         if (.not.havers) go to 125
         do 120 j=1,ncc
            t = c(i-1,j)
            c(i-1,j) = c(k,j)
  120    c(k,j) = t
  125    if (.not.wntv) go to 135
         do 130 j=1,nrv
            t = v(j,i-1)
            v(j,i-1) = v(j,k)
  130    v(j,k) = t
  135 continue
c                                  end of ordering algorithm.
  140 ier = 129
      if (fail) go to 9000
c                                  check for possible rank deficiency
      ier = 33
      t = 0.0
      if (d(1).ne.zero) t=d(n)/d(1)
      f=100.0+t
      if (f.eq.100.0) go to 9000
      ier = 0
      go to 9005
 9000 continue
      call uertst (ier,6hlsvdb )
 9005 return
      end
c   imsl routine name   - lsvg1
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - nucleus called only by imsl routine lsvdb
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine lsvg1  (a,b,cos,sin,sig)
c
c                                  specifications for arguments
      real               a,b,cos,sin,sig
c                                  specifications for local variables
      real               aa,bb
c                                  first executable statement
      if (abs(a).le.abs(b)) go to 5
      aa = abs(a+a)
      sig = aa*sqrt(0.25+(b/aa)**2)
      cos = a/sig
      sin = b/sig
      return
    5 if (b.eq.0.0) go to 10
      bb = abs(b+b)
      sig = bb*sqrt(0.25+(a/bb)**2)
      cos = a/sig
      sin = b/sig
      return
   10 sig = 0.0
      cos = 0.0
      sin = 1.0
      return
      end
c   imsl routine name   - lsvg2
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - nucleus called only by imsl routine lsvdb
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine lsvg2  (cos,sin,x,y)
c
c                                  specifications for arguments
      real               cos,sin,x,y
c                                  specifications for local variables
      real               xr
c                                  first executable statement
      xr=cos*x+sin*y
      y=-sin*x+cos*y
      x=xr
      return
      end
c   imsl routine name   - vhs12
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - real householder transformation -
c                           computation and applications.
c
c   usage               - call vhs12 (mode,lp,l1,m,u,incu,up,c,incc,
c                           icv,ncv)
c
c   arguments    mode   - option parameter. (input)
c                         if mode=1, the subroutine computes a
c                           householder transformation and if ncv.gt.0,
c                           multiplies it by the set of ncv vectors
c                           (each of length m) stored in c. for a
c                           given vector v of length m and two integer
c                           indices lp and l1 that satisfy
c                           1 .le. lp .lt. l1 .le. m, the subroutine
c                           defines an m by m householder
c                           transformation q which satifies qv=w where
c                           w(i)=v(i) for i.lt.lp
c                           w(lp)=-sig*sqrt(v(lp)**2+v(l1)**2+...
c                             +v(m)**2)
c                             sig=1  if v(lp).ge.0
c                             sig=-1 if v(lp).lt.0
c                           w(i)=v(i) for lp.lt.i.lt.l1
c                           w(i)=0    for i.ge.l1.
c                         if mode=2, the subroutine assumes that a
c                           householder transformation has already
c                           been defined by a previous call with
c                           mode=1, and if ncv.gt.0, multiplies it by
c                           the set of ncv vectors (each of length
c                           m) stored in c.
c                lp     - parameters that define the desired
c                l1         householder transformation. (input)
c                m          if the condition 1.le.lp.lt.l1.le.m is
c                           not satisfied, the subroutine returns to
c                           the calling program without performing
c                           any computations.
c                u      - vector of m elements. (input, and output if
c                           mode=1)
c                           the storage increment between elements
c                           of u is incu. (i.e., u(1+(j-1)*incu),
c                           j=1,...,m). if mode=1, the array v is
c                           defined as v(j)=u(1+(j-1)*incu),
c                           j=1,...,m.
c                         on output, u(1+(lp-1)*incu) is set to
c                           w(lp) (as defined above in the description
c                           of mode=1).
c                incu   - increment between elements of u. (input)
c                up     - scalar set to v(lp)-w(lp) to define the
c                           householder transformation q. (input if
c                           mode=2, output if mode=1)
c                c      - vector of ncv*m elements. (input/output)
c                           if ncv.le.0, c is not used.
c                           if ncv.gt.0, c contains ncv vectors
c                           of length m with increment incc between
c                           elements of vectors and increment icv
c                           between vectors. element i of vector j is
c                           defined as c(1+(i-1)*incc+(j-1)*icv),
c                           i=1,...,m and j=1,...,ncv.
c                         on output, c contains the set of ncv
c                           vectors resulting from multiplying
c                           the given vectors by q.
c                incc   - increment between elements of vectors
c                           in c. (input)
c                icv    - increment between vectors in c. (input)
c                ncv    - number of vectors stored in c. (input)
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks  1.  if u is a single subscripted array or the j-th column
c                of a matrix, then incu=1. if u is the i-th row of a
c                matrix then incu is the row dimension of the matrix
c                exactly as specified in the calling program.
c            2.  if c is a double subscripted matrix and the vectors
c                are the first ncv columns of c, then incc=1 and icv
c                is the row dimension of c exactly as specified in
c                the calling program. in this case c is replaced
c                by qc. if the vectors are successive rows of c
c                then incc is the row dimension of c exactly as
c                specified in the calling program and icv=1. in this
c                case c is replaced by cq.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine vhs12  (mode,lp,l1,m,u,incu,up,c,incc,icv,ncv)
c
c                                  specifications for arguments
      integer            mode,lp,l1,m,incu,incc,icv,ncv
      real               u(1),up,c(1)
c                                  specifications for local variables
      integer            ij,ilp,il1,im,incr,i2,i3,i4,j
      double precision   sm,b
      real               one,cl,clinv,sm1
c                                  first executable statement
      one = 1.
c
      if (0.ge.lp.or.lp.ge.l1.or.l1.gt.m) go to 9005
      ilp = (lp-1)*incu+1
      il1 = (l1-1)*incu+1
      im = (m-1)*incu+1
      cl = abs(u(ilp))
      if (mode.eq.2) go to 15
c                                  construct the transformation.
      do 5 ij=il1,im,incu
    5 cl = amax1(abs(u(ij)),cl)
      if (cl.le.0.0) go to 9005
      clinv = one/cl
      sm = (dble(u(ilp))*clinv)**2
      do 10 ij=il1,im,incu
   10 sm = sm+(dble(u(ij))*clinv)**2
c                                  convert dble. prec. sm to sngl.
c                                    prec. sm1
      sm1 = sm
      cl = cl*sqrt(sm1)
      if (u(ilp).gt.0.0) cl = -cl
      up = u(ilp)-cl
      u(ilp) = cl
      go to 20
c                                  apply the transformation
c                                    i+u*(u**t)/b to c.
   15 if (cl.le.0.0) go to 9005
   20 if (ncv.le.0) go to 9005
      b = dble(up)*u(ilp)
c                                  b must be nonpositive here. if b =
c                                    0., return.
      if (b.ge.0.0) go to 9005
      b = one/b
      i2 = 1-icv+incc*(lp-1)
      incr = incc*(l1-lp)
      do 35 j=1,ncv
         i2 = i2+icv
         i3 = i2+incr
         i4 = i3
         sm = c(i2)*dble(up)
         do 25 ij=il1,im,incu
            sm = sm+c(i3)*dble(u(ij))
            i3 = i3+incc
   25    continue
         if (sm.eq.0.0) go to 35
         sm = sm*b
         c(i2) = c(i2)+sm*dble(up)
         do 30 ij=il1,im,incu
            c(i4) = c(i4)+sm*dble(u(ij))
            i4 = i4+incc
   30    continue
   35 continue
 9005 return
      end
c   imsl routine name   - uerset
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - set message level for imsl routine uertst
c
c   usage               - call uerset (level,levold)
c
c   arguments    level  - new value for message level. (input)
c                           output from imsl routine uertst is
c                           controlled selectively as follows,
c                             level = 4 causes all messages to be
c                                       printed,
c                             level = 3 messages are printed if ier is
c                                       greater than 32,
c                             level = 2 messages are printed if ier is
c                                       greater than 64,
c                             level = 1 messages are printed if ier is
c                                       greater than 128,
c                             level = 0 all message printing is
c                                       suppressed.
c                levold - previous message level. (output)
c
c   precision/hardware  - single/all
c
c   reqd. imsl routines - uertst,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine uerset (level,levold)
c                                  specifications for arguments
      integer            level,levold
c                                  first executable statement
      levold = level
      call uertst (levold,6huerset)
      return
      end
c   imsl routine name   - uertst
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - december 1, 1980
c
c   purpose             - print a message reflecting an error condition
c
c   usage               - call uertst (ier,name)
c
c   arguments    ier    - error parameter. (input)
c                           ier = i+j where
c                             i = 128 implies terminal error,
c                             i =  64 implies warning with fix, and
c                             i =  32 implies warning.
c                             j = error code relevant to calling
c                                 routine.
c                name   - a six character literal string giving the
c                           name of the calling routine. (input)
c
c   precision/hardware  - single/all
c
c   reqd. imsl routines - ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks      the error message produced by uertst is written
c                onto the standard output unit. the output unit
c                number can be determined by calling ugetio as
c                follows..   call ugetio(1,nin,nout).
c                the output unit number can be changed by calling
c                ugetio as follows..
c                                nin = 0
c                                nout = new output unit number
c                                call ugetio(3,nin,nout)
c                see the ugetio document for more details.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine uertst (ier,name)
c                                  specifications for arguments
      integer            ier
      integer*2          name(3)
c                                  specifications for local variables
      integer*2          namset(3),nameq(3)
      data               namset/2hue,2hrs,2het/
      data               nameq/2h  ,2h  ,2h  /
c                                  first executable statement
      data               level/4/,ieqdf/0/,ieq/1h=/
      if (ier.gt.999) go to 25
      if (ier.lt.-32) go to 55
      if (ier.le.128) go to 5
      if (level.lt.1) go to 30
c                                  print terminal message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,35) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,35) ier,name
      go to 30
    5 if (ier.le.64) go to 10
      if (level.lt.2) go to 30
c                                  print warning with fix message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,40) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,40) ier,name
      go to 30
   10 if (ier.le.32) go to 15
c                                  print warning message
      if (level.lt.3) go to 30
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,45) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,45) ier,name
      go to 30
   15 continue
c                                  check for uerset call
      do 20 i=1,3
         if (name(i).ne.namset(i)) go to 25
   20 continue
      levold = level
      level = ier
      ier = levold
      if (level.lt.0) level = 4
      if (level.gt.4) level = 4
      go to 30
   25 continue
      if (level.lt.4) go to 30
c                                  print non-defined message
      call ugetio(1,nin,iounit)
      if (ieqdf.eq.1) write(iounit,50) ier,nameq,ieq,name
      if (ieqdf.eq.0) write(iounit,50) ier,name
   30 ieqdf = 0
      return
   35 format(19h *** terminal error,10x,7h(ier = ,i3,
     1       20h) from imsl routine ,3a2,a1,3a2)
   40 format(36h *** warning with fix error  (ier = ,i3,
     1       20h) from imsl routine ,3a2,a1,3a2)
   45 format(18h *** warning error,11x,7h(ier = ,i3,
     1       20h) from imsl routine ,3a2,a1,3a2)
   50 format(20h *** undefined error,9x,7h(ier = ,i5,
     1       20h) from imsl routine ,3a2,a1,3a2)
c                                  save p for p = r case
c                                    p is the page name
c                                    r is the routine name
   55 ieqdf = 1
      do 60 i=1,3
   60 nameq(i) = name(i)
   65 return
      end
c   imsl routine name   - ugetio
c
c-----------------------------------------------------------------------
c
c   computer            - prime/single
c
c   latest revision     - june 1, 1981
c
c   purpose             - to retrieve current values and to set new
c                           values for input and output unit
c                           identifiers.
c
c   usage               - call ugetio(iopt,nin,nout)
c
c   arguments    iopt   - option parameter. (input)
c                           if iopt=1, the current input and output
c                           unit identifier values are returned in nin
c                           and nout, respectively.
c                           if iopt=2, the internal value of nin is
c                           reset for subsequent use.
c                           if iopt=3, the internal value of nout is
c                           reset for subsequent use.
c                nin    - input unit identifier.
c                           output if iopt=1, input if iopt=2.
c                nout   - output unit identifier.
c                           output if iopt=1, input if iopt=3.
c
c   precision/hardware  - single/all
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c   remarks      each imsl routine that performs input and/or output
c                operations calls ugetio to obtain the current unit
c                identifier values. if ugetio is called with iopt=2 or
c                iopt=3, new unit identifier values are established.
c                subsequent input/output is performed on the new units.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - imsl warrants only that imsl testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c-----------------------------------------------------------------------
c
      subroutine ugetio(iopt,nin,nout)
c                                  specifications for arguments
      integer            iopt,nin,nout
c                                  specifications for local variables
      integer            nind,noutd
      data               nind/5/,noutd/6/
c                                  first executable statement
      if (iopt.eq.3) go to 10
      if (iopt.eq.2) go to 5
      if (iopt.ne.1) go to 9005
      nin = nind
      nout = noutd
      go to 9005
    5 nind = nin
      go to 9005
   10 noutd = nout
 9005 return
      end
