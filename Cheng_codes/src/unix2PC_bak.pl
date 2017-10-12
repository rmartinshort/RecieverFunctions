#! /usr/bin/perl
# perl scripts for convert files from UNIX format to PC format
# 09/10/08: I found sacsun2linux command does work in my machine !!! #can't work in my machine^^^^^^^^  so pity

#$DIR= "/home/lchen/c/D/RF/SRF/data/bhw_eqz_UNIX";
#$DIRNEW= "/home/lchen/c/D/RF/SRF/data/bhw_eqz";
$DIR= "/home/weizg/test";
$DIRNEW= "/home/weizg/hbb_temppc";

(-d $DIRNEW) || mkdir($DIRNEW,0744);
chdir $DIR;

@DIRS = `ls `;

for($j = 0; $j < @DIRS; $j++){
 $dir = $DIRS[$j]; chomp($dir);

# if( -d $dir){
 if( -d $dir && $dir =~ /2/){

  chdir $dir;

  mkdir("$DIRNEW/$dir", 0744);

  print stdout "changing to directory $dir\n";

  @STATION = `ls *.* `;
 for($st = 0; $st < @STATION; $st++ ){
   $SACFILE = $STATION[$st]; chomp($SACFILE);
   print STDOUT "the sacfile is $SACFILE\n";	
  `cp $SACFILE $DIRNEW/$dir/$SACFILE`;
  `sacsun2linux $DIRNEW/$dir/$SACFILE`;
  }
  chdir("..");
 }
}


