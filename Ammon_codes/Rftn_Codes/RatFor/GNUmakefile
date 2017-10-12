#
#DEBUG := yes
# pd ratfor77 (oz)
#
# if F77 is defined, the output
# of ratfor is Fortran 77.
#
#	On sun4,		use S_CHAR="char"
#	On RS6000,		use S_CHAR="signed char"
#	On DEC3100,	maybe	use S_CHAR="signed char"
#	On CRAY,		use S_CHAR="char"
#	On GNU,		        use S_CHAR="char"
#
# Default definition of the makefile include files
# If you have kept the whole distribution together you won't need to set
# this yourself. If you have just taken a few directories you should set
# the envionment variable "MAKEINC" to point at the config directory.
#
ifndef MAKEINC
    MAKEINC := ./../../config
endif
#

include $(MAKEINC)/Makefile.SEP.defs

# we don't want the standard SEP flags so override them
ifeq ($(GNU),yes)
CFLAGS := -c -DF77 -DS_CHAR=$(SIGNED_CHAR) -DGNU
else
CFLAGS := -c -DF77 -DS_CHAR=$(SIGNED_CHAR)
endif

all: $(SEPBINDIR)/ratfor77
	@echo "making all in ratfor77 : done"
#
deinstall:
	$(RM) $(SEPBINDIR)/ratfor77
	@echo "deinstall in ratfor77 : done"
#
$(SEPBINDIR)/ratfor77:	$(addprefix $(MACHINETYPE)/, rat4.o lookup.o getopt.o )
	 $(link.c)
	$(INSTALL) $(INSTALLBIN) $(MACHINETYPE)/a.out $@

include $(MAKEINC)/Makefile.rules

