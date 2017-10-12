        subroutine  fork(lx,cx,signi)
c  **********************************************************************
c  *  FAST FOURIER TRANSFORMATION                                       *
c  *  input output cx(lx) complex                                       *
c  *  lx is the interger power to 2                                     *
c  *  signi -1 forward fft                                              *
c  *        +1 inverse fft                                              *
c  *                       lx                                           *
c  *  cx(k) = sqrt(1./lx) sum (cx(j)*exp(2*pi*signi*i*(j-1)*(k-1)/lx))  *
c  *                      j=1                                           *
c  *       for k=1,2,...,lx                                             *
c  **********************************************************************
        complex cx(lx),carg,cexp,cw,ctemp
        integer signi
        j=1
        sc=1.
        if(signi.eq.1) sc=1./float(lx)
c        sc=sqrt(1./lx)
        do 30 i=1,lx
        if(i.gt.j)  go to 10
        ctemp=cx(j)*sc
        cx(j)=cx(i)*sc
        cx(i)=ctemp
  10    m=lx/2
  20    if(j.le.m) go to 30
        j=j-m
        m=m/2
        if(m.ge.1) go to 20
  30    j=j+m
        l=1
  40    istep=2*l
        do 50 m=1,l
        carg=(0.,1.)*(3.14159265*signi*(m-1))/l
        cw=cexp(carg)
        do 50 i=m,lx,istep
        ctemp=cw*cx(i+l)
        cx(i+l)=cx(i)-ctemp
  50    cx(i)=cx(i)+ctemp
        l=istep
        if(l.lt.lx) go to  40
        return
        end
