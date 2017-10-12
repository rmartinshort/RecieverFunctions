#! /bin/perl
# perl scripts for filter eqr and eqt files      ## can't filter 0.4 to 0.3  it will produce sharp-angle , but 1 is ok  smooth-angle 2008/12/23
# Notice if you need to do decimation

$DIR = "/home/weizg/testlinux";
 $DIRTMP= "/home/weizg/f3"; 

@fmin = (0.03,0.03,0.03); @fmax = (0.3,0.5,1);  #   2008/11/17
$i=0;
$appendix="eq*";
$SUBDIR="202";
(-d $DIRTMP) || mkdir($DIRTMP,0744); 


if( -d $DIR){
  chdir $DIR;

  @DIRS = `ls `;

  for($j = 0; $j < @DIRS; $j++){
    $dir = $DIRS[$j]; chomp($dir);

    if( -d $dir){
#    if( -d $dir && $dir =~ /$SUBDIR/){

      print stdout "changing to directory $dir\n";
      chdir $dir;
      (-d "$DIRTMP/$dir") || mkdir("$DIRTMP/$dir", 0755);
    
      @STATION = `ls *.$appendix `;
      for($st = 0; $st < @STATION; $st++ ){
        $SACFILE = $STATION[$st]; chomp($SACFILE);

       print STDOUT "the sacfile is $SACFILE\n";	


    $fmin=$fmin[$i];
    $fmax=$fmax[$i];
    open( SAC, "|sac2000");       
    print SAC "r $SACFILE\n";  # shu chu dao SAC  2008/11/17
#    print SAC "rmean\n";
#    print SAC "rtrend\n";
#    print SAC "taper\n";
    print SAC "bp bu n 2 p 2 co $fmin $fmax\n";	  #bandpass   sac shouce  p59
#    print SAC "trans from polezero subtype $RESPONSE to none\n";
#    print SAC "bp co $fmin $fmax n 2 p 2\n";
#    print SAC "decimate 4\n";   # sac cankao   p82   meigesigeshujuquyigezhi  byme 2009/5/9
    print SAC "w $DIRTMP/$dir/$SACFILE\n";
    print SAC "quit\n";
    close(SAC);

    print STDOUT "$SACFILE was filtered and copied\n";
      }
      chdir("..");
    }
  }
}
