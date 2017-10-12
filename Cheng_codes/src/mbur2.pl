#! /usr/bin/perl
# perl scripts for calculating the receiver function with 
# time-domain deconvolution

my($onesta)=@ARGV;
#$DIR= "/home/lchen/c/40GDisk_bak/disk-1/rf/data_a2my/bh1";
#$DIR= "/home/lchen/RF/data/data_a2my/taiwan";/add a denotation "#"(08/10/7)/
$DIR="/home/weizg/data/ncisp521";
chdir($DIR);

@DIRS = `ls `;
for($j = 0; $j < @DIRS; $j++){ 
 $dir = $DIRS[$j]; chomp($dir);
 if($onesta ne ""){$dir=$onesta;chomp($dir);}
 $num = $dir;
 if( -d $dir ){
 #if( -d $dir && $dir =~ /^c$/){print"hello/n";
# if( -d $dir && $dir =~ /02[\d]/){hbb_temppc
# if(-d $dir && ($dir =~ /^0[89]/ || $dir =~ /10[\d]/)){
# if( -d $dir && $num < 2){

  chdir $dir;
  print stdout "*********************************************\n";
  print stdout "*DIRECTORY:$DIR/$dir\n";
  print stdout "*********************************************\n";

  !(-f "record") || (unlink "record");
  !(-f "rstation") || (unlink "rstation");print"hello/n";
  `ls *.bhe > record`;    # > record: input to record (08/10/7)#print"hello/n";
  `awk \-F\. \'\{print \$1\}\' record \> rstation`;#input sth into rstation 08/10/7#
#(08/10/7 change this conmand line into note line ) system("/home/lchen/program/program/decon/xww/burgays1_test");
system("/home/weizg/plfc/RF_calculation/burgays1_test");

  chdir("..");
  if($dir eq $onesta){die("The directory $dir has finished!\n");}
 }
}

