#RMS 2017
#Some of this is not used 
#infact, its likely that only the burgays code will be useful. This is the only one that we 
#compile at the moment

#Here, libs points to the SAC library. Note that for these codes to work the data needs first
#to be converted into RTZ 


LIBS = -L/home/rmartinshort/Software/sac/sac/lib -lm -lsac ./libsac.a  ./libsacio.a  ./sacio.a
F77  = gfortran -ffixed-line-length-0
BIN  = ../bin

ALL: burgays

burgays =  ask.o asktxt.o blank.o \
    minmax.o  burgays.o\
    yesno.o zero.o dfftr.o npowr2.o \
    fft.o locast.o rotate.o getdata.o filter.o jacobi.o canc.o cross.o \
    convolution.o lowpass.o gaussfilter.o isvd.o csoln.o
burgays: $(burgays)
	$(F77) -g $(burgays) -o $(BIN)/burgays $(LIBS)

burgays1_test =  ask.o asktxt.o blank.o \
    minmax.o  burgays1_test.o\
    yesno.o zero.o dfftr.o npowr2.o \
    fft.o locast.o rotate.o getdata.o filter.o jacobi.o canc.o cross.o \
    convolution.o lowpass.o gaussfilter.o isvd.o csoln.o
burgays1_test: $(burgays1_test)
	$(F77) -g $(burgays1_test) -o $(BIN)burgays1_test $(LIBS)

burgays_zrt1 =  ask.o asktxt.o blank.o \
    minmax.o  burgays_zrt1.o\
    yesno.o zero.o dfftr.o npowr2.o \
    fft.o locast.o rotate.o getdata.o filter.o jacobi.o canc.o cross.o \
    convolution.o lowpass.o gaussfilter.o isvd.o csoln.o
burgays_zrt1: $(burgays_zrt1)
	$(F77) -g $(burgays_zrt1) -o $(BIN)burgays_zrt1 $(LIBS)

burgays_s =  ask.o asktxt.o blank.o \
    minmax.o  burgays_s.o\
    yesno.o zero.o dfftr.o npowr2.o \
    fft.o locast.o rotate.o getdata.o filter.o jacobi.o canc.o cross.o \
    convolution.o lowpass.o gaussfilter.o isvd.o csoln.o
burgays_s: $(burgays_s)
	$(F77) -g $(burgays_s) -o $(BIN)burgays_s $(LIBS)

ask.o: \
	ask.f
	$(F77) -c ask.f
asktxt.o: \
	asktxt.f
	$(F77) -c asktxt.f
blank.o: \
	blank.f
	$(F77) -c blank.f
minmax.o: \
	minmax.f
	$(F77) -c minmax.f
burgays.o: \
	burgays.f
	$(F77) -c burgays.f
burgays1_test.o: \
	burgays1_test.f
	$(F77) -c burgays1_test.f
burgays_zrt1.o: \
	burgays_zrt1.f
	$(F77) -c burgays_zrt1.f
burgays_s.o: \
	burgays_s.f
	$(F77) -c burgays_s.f
yesno.o: \
	yesno.f
	$(F77) -c yesno.f
zero.o: \
	zero.f
	$(F77) -c zero.f
dfftr.o: \
	dfftr.f
	$(F77) -c dfftr.f
npowr2.o: \
	npowr2.f
	$(F77) -c npowr2.f
fft.o: \
	fft.f
	$(F77) -c fft.f
locast.o: \
	locast.f
	$(F77) -c locast.f
rotate.o: \
	rotate.f
	$(F77) -c rotate.f
getdata.o: \
	getdata.f
	$(F77) -c getdata.f
filter.o: \
	filter.f
	$(F77) -c filter.f
jacobi.o: \
	jacobi.f
	$(F77) -c jacobi.f
canc.o: \
	canc.f
	$(F77) -c canc.f
cross.o: \
	cross.f
	$(F77) -c cross.f
convolution.o: \
	convolution.f
	$(F77) -c convolution.f
lowpass.o: \
	lowpass.f
	$(F77) -c lowpass.f
gaussfilter.o: \
	gaussfilter.f
	$(F77) -c gaussfilter.f
isvd.o: \
	isvd.f
	$(F77) -c isvd.f
csoln.o: \
	csoln.f
	$(F77) -c csoln.f

clean:
	rm -f *.o
