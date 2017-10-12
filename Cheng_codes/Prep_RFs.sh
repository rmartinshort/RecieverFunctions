#!/bin/bash
#RMS Feb 2017
#Following Cheng's workflow in preparing data for RF analysis

#Enter directory structure created by obspyDMT, make a new directory for the RFs, copy the
#SAC files here, filter at the desired passband.

freqmin=2
freqmax=4

for event in $( ls -d 20* ); do 

echo In $event
	
cd $event
cd BH_VEL
rm -r rf
mkdir -p rf

cp *BHZ rf
cp *BHR rf
cp *BHT rf

cd rf

#-------------------------------------------------------------------------------------
#look at each sacfile, filter in the desired passband and cut around the P arrival

for sacfile in $( ls *BH* ); do

echo $sacfile
parrival=`saclst USER1 f $sacfile | awk '{print $2}'`
parr_lower=`echo $parrival | awk -F" " '{print ($1-60)}'`
parr_upper=`echo $parrival | awk -F" " '{print ($1+80)}'`


#--------------------------------------------------------------------------------------
#In this case, we're limiting the upper passband to 1/freqmin seconds (typically 5 or 10?) 

sac<<EOF
cut $parr_lower $parr_upper
read $sacfile
taper w 0.1
rmean
rtrend
bp co 0.01 0.1 p $freqmin n $freqmax
write $sacfile.filt
quit
EOF

rm $sacfile

done

#--------------------------------------------------

#Creates a list of unique station names, followed by the character 'n'
#This is needed as input to the program 'burgays'
#uses the uniq program, which outputs all lines exactly once (avoids repeats)
ls *BHZ.filt | awk -F. '{print $3}' | uniq > rstation
echo n >> rstation

cd ../../../

done 
