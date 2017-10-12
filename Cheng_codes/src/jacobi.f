       subroutine jacobi(n,ev,a,d,v,b,z,rt)
       dimension a(n,n),d(n),v(n,n),b(n),z(n)
       logical ev
       integer p,q,p01,q01,p11,q11,rt
       n11=n-1
       if(.not.ev) go to 104
       do 103 p=1,n
       do 103 q=1,n
       if((p-q)) 102,101,102
 101   v(p,q)=1.0
       go to 103
 102   v(p,q)=0.0
 103   continue
 104   do 105 p=1,n
       b(p)=a(p,p)
       d(p)=a(p,p)
       z(p)=0.0
 105   continue
       rt=0.0
       do 122 i=1,50
       sm=0.0
       do 106 p=1,n11
       p01=p+1
       do 106 q=p01,n
       sm=sm+abs(a(p,q))
 106   continue
       if(sm.eq.0.) return
       if(i-4) 107,108,108
 107   tres=0.2*sm/(n*n)
       go to 109
 108   tres=0.0
 109   do 120 p=1,n11
       p01=p+1
       do 120 q=p01,n
       g=100*abs(a(p,q))
       if(.not.(i.gt.4.and.abs(d(p))+g.eq.abs(d(p)).and.abs(d(q))
     *  +g.eq.abs(d(q)))) go to 110
       a(p,q)=0
       go to 120
 110   if(abs(a(p,q)).le.tres) go to 120
       h=d(q)-d(p)
       if(abs(h)+g.eq.abs(h)) go to 111
       thet=0.5*h/a(p,q)
       t1=abs(thet)
       t2=sqrt(t1)*sqrt(1./t1+t1)
       t=1./(t1+t2)
       if(thet.lt.0.) t=-1.*t
       go to 112
 111   t=a(p,q)/h
 112   c=1./sqrt(1.+t**2)
       s=t*c
       tau=s/(1+c)
       h=t*a(p,q)
       z(p)=z(p)-h
       z(q)=z(q)+h
       d(p)=d(p)-h
       d(q)=d(q)+h
       a(p,q)=0.
       p11=p-1
       if(p11.lt.1) go to 214
       do 113 j=1,p11
       g=a(j,p)
       h=a(j,q)
       a(j,p)=g-s*(h+g*tau)
       a(j,q)=h+s*(g-h*tau)
 113   continue
 214   q11=q-1
       if(q11.lt.p01) go to 115
       do 114 j=p01,q11
       g=a(p,j)
       h=a(j,q)
       a(p,j)=g-s*(h+g*tau)
       a(j,q)=h+s*(g-h*tau)
 114   continue
 115   q01=q+1
       if(q01.gt.n) go to 119
       do 116 j=q01,n
       g=a(p,j)
       h=a(q,j)
       a(p,j)=g-s*(h+g*tau)
       a(q,j)=h+s*(g-h*tau)
 116   continue
 119   if(.not.ev) go to 118
       do 117 j=1,n
       g=v(j,p)
       h=v(j,q)
       v(j,p)=g-s*(h+g*tau)
       v(j,q)=h+s*(g-h*tau)
 117   continue
 118   rt=rt+1
 120   continue
       do 121 p=1,n
       d(p)=b(p)+z(p)
       b(p)=d(p)
       z(p)=0.
 121   continue
 122   continue
       return 
       end
        
