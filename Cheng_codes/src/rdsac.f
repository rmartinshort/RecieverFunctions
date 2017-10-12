c     read a sac file data
      subroutine rdsac(file1,x)
      dimension x(1)
      include 'hdr.inc'
      character*32 file1
      inunit1=8
      open(inunit1,file=file1,access='direct',recl=158*4)
      read(inunit1,rec=1) (rheader(i),i=1,70),(iheader(i),i=1,35),
     *  (lheader(i),i=1,5),kstnm,kevnm,(cheader(i),i=1,21)
      close(inunit1)
      k=(158+iheader(10))*4
      open(inunit1,file=file1,access='direct',recl=k)
      read(inunit1,rec=1) (idummy,i=1,158),(x(i),i=1,iheader(10))
      close(inunit1)
      return
      end
