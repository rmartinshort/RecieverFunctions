      subroutine whdr
      include 'hdr.inc'
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c
c **************************************************************************
c
c    call to whdr is:
c                      call whdr
c
c
c ******************************************************************************
c
c
c     write a sac header       
C
c
      do 10 i=1,70
      rheader(i)=-12345.
   10 continue
      do 20 i=1,35
      iheader(i)=-12345
   20 continue
      do 30 i=1,21
      cheader(i)='-12345'
   30 continue
      do 40 i=1,5
      lheader(i)=.false.
   40 continue
      kstnm='-12345'
      kevnm='-12345'
      inunit1=8
      iheader(7)=6
      iheader(16)=1
      lheader(1)=.true.
      lheader(3)=.true.
      lheader(4)=.true.
      length=158*4
      return
      end
