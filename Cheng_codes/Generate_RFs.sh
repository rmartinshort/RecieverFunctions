#!/bin/bash
#RMS Feb 2017
#Following Cheng's workflow in preparing data for RF analysis 
#This appears to generate the reciever functions using the 'burgays' code

for event in $( ls -d 20* ); do

    cd $event/BH_VEL/rf/

    rm *.eqr *.eqz *.eqt

    for sacfile in $( ls *..BHR.filt ); do
	filename1=`echo $sacfile |  awk -F. '{print $3}'`
	filename2=$filename1.r
	mv $sacfile $filename2
    done

    for sacfile in $( ls *..BHZ.filt ); do
	filename1=`echo $sacfile | awk -F. '{print $3}'`
        filename2=$filename1.z
	mv $sacfile $filename2
    done

    for sacfile in $( ls *..BHT.filt ); do
        filename1=`echo $sacfile | awk -F. '{print $3}'`
	filename2=$filename1.t
	mv $sacfile $filename2
    done

sac<<EOF
r *.z
decimate 2
w over
quit
EOF

#Run the burgays code, which will generate the receiver functions
/home/rmartinshort/Documents/Berkeley/Receiver_functions/simple_stack/RF_calc/burgays rstation

cd ../../../

done

