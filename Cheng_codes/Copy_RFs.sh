#!/bin/bash
#RMS Feb 2017
#Following Cheng's RF workflow

#Copy the reciever functions produced in each of the event directories to a new folder, where they can undergo further processing


mkdir -p RFdata

for event in $( ls -d 20*/ ); do

    cd $event/BH_VEL/rf/

    for sacfile in $( ls *eqr ); do

	npts=`saclst NPTS f $sacfile | awk '{print $2}'`
	eve=`echo $event | awk -F/ '{print $1}'`
	filename2=$eve.$sacfile

	if [ $npts -le 10001 ]; then
sac<<EOF
r $sacfile
write over
quit
EOF
	    cp $sacfile ../../../RFdata/$filename2
	fi
    done

    cd ../../../

    cd RFdata
    ls *.eqr > data
    cd ../

done
