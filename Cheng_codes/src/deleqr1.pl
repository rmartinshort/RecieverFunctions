#! /bin/perl
# perl scripts for deleting eqrfile. 
#Usage:perl $0 ($num $sta)



#$DIR = "/work1/common/home/chl/RF/Japan_pPf05";
#$DIR = "/work1/common/home/chl/RF/hbbf05";
#$DIR = "/work1/common/home/chl/SRF/hbbf03_dt01";
#$DIR = "/work1/common/home/chl/SRF/SDQ1f03_dt01";
#$DIR = "/work1/common/home/chl/SRF/bhwf03_dt01";
 $DIR="/home/weizg/data/test05";
#$DIR = "/work1/common/home/chl/SRF/Aif03_dt01";
#$DIR = "/work1/common/home/chl/SRF/wnccf03_dt01";
$dirbad="bad1";

(-d $DIR) || die("The directory $DIR doesnot exist!\n");
my($rec,$num,$head)=@ARGV;
if($num eq ""){$num=9;}
chdir($DIR);
@DIRS = `ls `;
for($j = 0; $j < @DIRS; $j++){
 $dir = $DIRS[$j]; chomp($dir);
 if($rec ne ""){$dir=$rec;}
 if( -d $dir){

  chdir($dir);
  print stdout "changing to directory $dir\n";
  $DIRNOW = `pwd`; chomp($DIRNOW);
  $newdir = $DIRNOW."/".$dirbad;    #make a bad-dir in /2/
  (-d $newdir) || mkdir($newdir,0755);
  # chdir(bad); # add by myself for choosing available data from bad the first time

  @STATION = `ls *.eqr`;  # @STATION = `ls *.eqz `; 

#sorting by the $head.
  if($head eq ""){$head="az";}
  print stdout "Begin sorting\n";
  for($si=0; $si<@STATION; $si++){ 
    $hdvalue[$si]=`saclst $head f $STATION[$si]`;  chomp($STATION[$si]);# `saclst $head -h <$STATION[$si]`;
   ($temp1,$hdvalue[$si]) = split(' ',$hdvalue[$si]); chomp($stla); 
   print stdout "$si:$STATION[$si]::$hdvalue[$si]\n";
  }
  $sn=@STATION;
  &pcbub(@hdvalue,@STATION,$sn); # &pcbub(\@hdvalue,\@STATION,$sn);yexing?
  print stdout "After sorting\n";
  for($sl=0; $sl<@STATION; $sl++){ 
    print stdout "$sl:$STATION[$sl]::$hdvalue[$sl]\n";
  }
#sorted by the az. 

  for($st = 0; $st < @STATION; $st=$st+$num){
    $jst=$st;
        for($jj=0; $jj < $num; $jj++,$jst++){
          $sacarray[$jj]=$STATION[$jst];chomp($sacarray[$jj]);
        }
	open(SAC,"|sac2000");
#print stdout "\n**********r @sacarray****\n\n";
	print SAC "cut -10 100\n";
	print SAC "r @sacarray\n";
        print SAC "fileid format equals type list kstnm az\n";
#        print SAC "gtext size large\n";
        print SAC "qdp off\n";
        print SAC "ppk\n";
        print SAC "ppk\n";
	print SAC "quit\n";
	close(SAC);
        print stdout "\nFile list:\n";
        for($ii=1; $ii<=$num; $ii++){
           print stdout "$ii:$sacarray[$ii-1]\n";}
        print stdout "Which one you want to delete(1,2,3...)\n";
        $flag = <stdin>;
        if(($flag =~ /[^(\d&\s)]/)||($flag =~ /^S/))
            {die("Program terminated(by user)!\n");}
	for($jj=1; $jj<=$num; $jj++){
          if($flag =~ /$jj/){
	     ($sacf,$sactail)=split(/\./,$sacarray[$jj-1]);
	     if($sacf eq ""){last;}
             print stdout "Now deleting $sacf*...\n";
             `mv $sacf* $newdir`;}
          if($flag =~ /^quit/){die("Program terminated(quit)!\n");}
        }
   }
  chdir("..");
  if($dir eq $rec){die("The directory $rec has finished!\n");}
 }
}

sub pcbub{
  my($ps,$names,$ns)=@_;
# $ns=@$ps;
  $ks=0; $ms=$ns-1;
  while ($ks<$ms)
    { $js=$ms-1; $ms=0;
      for ($is=$ks; $is<=$js; $is++){
        if ($$ps[$is]>$$ps[$is+1])
          { $ds=$$ps[$is]; $$ps[$is]=$$ps[$is+1]; $$ps[$is+1]=$ds; $ms=$is;
            $tf1=$$names[$is]; $$names[$is]=$$names[$is+1]; $$names[$is+1]=$tf1;}}
      $js=$ks+1; $ks=0;
      for ($is=$ms; $is>=$js; $is--){
        if ($$ps[$is-1]>$$ps[$is])
          { $ds=$$ps[$is]; $$ps[$is]=$$ps[$is-1]; $$ps[$is-1]=$ds; $ks=$is;
            $tf2=$$names[$is]; $$names[$is]=$$names[$is-1]; $$names[$is-1]=$tf2;}}
    }
  return(@$names);
}
