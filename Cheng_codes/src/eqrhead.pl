#! /usr/bin/perl
# perl scripts for add head to eqr files. 

$DIR= "/home/weizg/testlinux";

my($onesta)=@ARGV;
unless(open(ppofile,"nmo.dat")){
die("cannot open input file nmo.dat!\n");
}
@ppo = <ppofile>;
for($i=0; $i < @ppo; $i++){
 chomp($ppo[$i]);
 ($gc[$i],$pp[$i])=split(/,/,$ppo[$i]);
  $pp[$i]=$pp[$i]/111.195*1000000;
  ($pp[$i],$junk)=split(/\./,$pp[$i]);
  $pp[$i]=$pp[$i]/1000000;
#print stdout "pppppppppppppppp$pp[$i]\n";
}
close(ppofile);
$num = 1;


chdir $DIR;
@DIRS = `ls `;
for($j = 0; $j < @DIRS; $j++){
 $dir = $DIRS[$j]; chomp($dir);

 if($onesta ne ""){$dir=$onesta;chomp($dir);}
# if( -d $dir && $dir =~ /[\d]/){
# if( -d $dir && $dir =~ /029/){
 if( -d $dir ){

  chdir $dir;
  
  print stdout "changing to directory $dir\n";

  @STATION = `ls *.eqr`;

  for($st = 0; $st < @STATION; $st++ ){
    $SACFILE = $STATION[$st]; chomp($SACFILE);
    print STDOUT "the sacfile is $SACFILE\n";	    
    
    $temp = `saclst stla f $SACFILE `; chomp($temp);
    ($temp1,$stla) = split(' ',$temp); chomp($stla);
    $temp = `saclst stlo f $SACFILE `; chomp($temp);
    ($temp1,$stlo) = split(' ',$temp); chomp($stlo);
    $temp = `saclst stel f $SACFILE `; chomp($temp);
    ($temp1,$stel) = split(' ',$temp); chomp($stel);
    $temp = `saclst stdp f $SACFILE `; chomp($temp);
    ($temp1,$stdp) = split(' ',$temp); chomp($stdp);
    $temp = `saclst evla f $SACFILE `; chomp($temp);
    ($temp1,$evla) = split(' ',$temp); chomp($evla);
    $temp = `saclst evlo f $SACFILE `; chomp($temp);
    ($temp1,$evlo) = split(' ',$temp); chomp($evlo);
    $temp = `saclst evdp f $SACFILE `; chomp($temp);
    ($temp1,$evdp) = split(' ',$temp); chomp($evdp);
    $temp = `saclst kstnm f $SACFILE `; chomp($temp);
    ($temp1,$kstnm) = split(' ',$temp); chomp($kstnm);
    $temp = `saclst gcarc f $SACFILE `; chomp($temp);
    ($temp1,$gcarc) = split(' ',$temp); chomp($gcarc);
    ($temp1,$temp2) = split(/\./,$SACFILE);

    $eqrfile = join(/\./,$temp1,".eqr"); chomp($eqrfile); 
    $eqtfile = join(/\./,$temp1,".eqt"); chomp($eqtfile); 
    for($m=0; $m < @gc; $m++){
     if(abs($gc[$m]-$gcarc)<=0.25){$user0=$pp[$m];last;}
     } 

    ($amarker,$p)=&timeref($evdp,$gcarc,$num); 
    $p=$p/111.195*1000000;
    ($p,$junk)=split(/\./,$p);
    $p=$p/1000000;
    for($m=1; $m < @pp; $m++){
      if($p>=$pp[$m]){
        $user1=$pp[$m];
        $midd=($pp[$m-1]+$pp[$m])*0.5;
        if($p > $midd){$user1=$pp[$m-1];}
        last;
      }
    }

    open(SAC, "|sac2000");
    print SAC "r $eqrfile $eqtfile\n";
    print SAC "ch stla $stla\n";	 
    print SAC "ch stlo $stlo\n";	 
    print SAC "ch stel $stel\n";	 
    print SAC "ch stdp $stdp\n";	 
    print SAC "ch evla $evla\n";	 
    print SAC "ch evlo $evlo\n";	 
    print SAC "ch evdp $evdp\n";	 
    print SAC "ch kstnm $kstnm\n";	 
    print SAC "ch user0 $user0\n";	 
    print SAC "ch user1 $user1\n";	 
#print stdout "eqrfile########################:$eqrfile\n";
    print stdout "user0#######################:$user0\n";
    print SAC "w over\n";
    print SAC "quit\n";
    close(SAC);
    

    print STDOUT "The head of $eqrfile $eqtfile was added!\n";
  }
 chdir("..");
 if($onesta eq $dir){die("The directory has finished!\n");}
}
}


sub timeref{
    my ($sdepth,$gcarc,$num) = @_;
        ($tref,$p) = &ttimes($sdepth,$gcarc,$num);
    return($tref,$p);
}


sub ttimes{
  my ($depth,$gcarc,$num)=@_;
  $tmp_file = "/home/weizg/c";
  my @ll;
  my $l;
  my $t;
  my @c;
  open(TT,"|ttimes > $tmp_file") || die "Cannot start ttimes!\n";
  print TT "all\n\n$depth\n$gcarc\n-1\n-1\n";  close(TT);
  open(TT,"< $tmp_file");  @ll = <TT>;  close(TT);
  unlink $tmp_file;
  $t = -1;
#  print stdout "num = $num\n";
  while ($l = shift @ll) {
    if($num <= 1) {
    if($l =~ /^Source depth \(km\):/){
       shift @ll;       shift @ll;
       @c = split(' ',shift @ll);
       return($c[3],$c[6]);
       print stdout "ERROR in computing the delay time!\n";
    }
    }
    elsif($num == 2) {
    if($l =~ /^Source depth \(km\):/){
       shift @ll;       shift @ll;
       while ($l = @ll) {
         @c = split(' ',shift @ll);
#         print stdout "$c[0], $c[1], $c[2]\n";
         if($c[1] =~ /^PP$/) {
           print stdout "$num, $c[0], $c[1], $c[2], $c[5]\n";
           return($c[2],$c[5]);
           print stdout "ERROR in computing the delay time!\n";
         }
       }
    }
    }
    elsif($num == 3) {
    if($l =~ /^Source depth \(km\):/){
       shift @ll;       shift @ll;
       while ($l = @ll) {
         @c = split(' ',shift @ll);
#         print stdout "$c[0], $c[1], $c[2]\n";
         if($c[1] =~ /^pP$/) {
           print stdout "$num, $c[0], $c[1], $c[2], $c[5]\n";
           return($c[2],$c[5]);
           print stdout "ERROR in computing the delay time!\n";
         }
       }
    }
    }

  }
}
