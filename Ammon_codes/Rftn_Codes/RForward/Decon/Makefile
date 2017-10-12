.SUFFIXES: .r
.r.f:
	ratfor -C $< > $@
FFLAGS = -g
#To optimize, use the next FFLAGS
#FFLAGS = -O4
TJOLIB = ../Subs/subs.a
BINDIR = ../../bin

pwaveqn: pwaveqn.f	
	f77 $(FFLAGS) -o $(BINDIR)/pwaveqn pwaveqn.f $(TJOLIB) $(SACLIB)

