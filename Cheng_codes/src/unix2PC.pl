#! /usr/bin/perl
# perl scripts for convert files from UNIX format to PC format    this can work in my machine!weizigen  2008/11/27  
# 09/10/08: I found sacsun2linux command does work in my machine !!!
#           I try sacswap        

$DIR= "/home/weizg/test/";
$DIRNEW= "/home/weizg/testlinux";
(-d $DIRNEW) || mkdir($DIRNEW,0744);
chdir $DIR;

@DIRS = `ls`;
print "@DIRS hello\n";

for($j = 0; $j < @DIRS; $j++){
 $dir = $DIRS[$j]; chomp($dir);

 if( -d $dir){
# if( -d $dir && $dir =~ /21/){

   chdir $dir;

  mkdir("$DIRNEW/$dir", 0744);

  print stdout "changing to directory $dir\n";

 @STATION = `ls `;#`ls *.SAC`;
 for($st = 0; $st < @STATION; $st++ ){
    $SACFILE = $STATION[$st]; chomp($SACFILE);
    $SACFILE1 = $SACFILE.'.swap';
    print STDOUT "the sacfile is $SACFILE\n";	
   `cp $SACFILE $DIRNEW/$dir/$SACFILE`;
   `/home/weizg/plfc/RF_calculation/utils/sacswap $DIRNEW/$dir/$SACFILE`;
   `mv $DIRNEW/$dir/$SACFILE1 $DIRNEW/$dir/$SACFILE`;
  }
  chdir("..");
 }
}


