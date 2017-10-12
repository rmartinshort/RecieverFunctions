FFLAGS = -g
#To optimize, use the next FFLAGS
#FFLAGS = -O4
SACLIBRARY = $(SACLIB)
TJOLIB = ../Subs/subs.a
BINDIR = ../../bin

OBJ = respknt.o kntfun.o rcvrfn.o ifmat.o abm.o
respknt: $(OBJ) kennet.inc
	f77 $(FFLAGS) -o $(BINDIR)/respknt  $(OBJ) $(TJOLIB) $(SACLIBRARY)
respknt.o: respknt.f kennet.inc
rcvrfn.o: rcvrfn.f kennet.inc
ifmat.o:  ifmat.f kennet.inc
kntfun.o: kntfun.f 
abm.o: abm.f

clean:
	/bin/rm *.o
