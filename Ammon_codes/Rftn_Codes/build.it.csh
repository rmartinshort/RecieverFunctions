#!/bin/csh
#
#  EDIT THE NEXT LINE
#
setenv SACLIB /Volumes/hd1/Users/cammon/Programs/CJA_Sacio/bin/isacio.a
#
echo ' '
echo 'Making the bin directory'
echo ' '
mkdir ./bin
#
cd RForward
#
echo ' '
echo 'Building the library'
echo ' '
cd Subs
make
make clean
#
echo ' '
echo 'Compiling the deconvolution program'
echo ' '
cd ../Decon
make
#
echo ' '
echo 'Compiling the forward-modeling program'
echo ' '
cd ../RespKennett
make
make clean
#
echo ' '
echo 'Compiling the common utilities'
echo ' '
cd ../Utilities
make mostcommon
#
echo ' '
echo 'Compiling the inversion programs'
echo ' '
cd ../../RInversion
make all
make clean
#
echo ' '
echo 'Executables are in the directory ./bin'
echo ' '
