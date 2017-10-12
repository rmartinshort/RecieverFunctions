########## 2008/11/21 weizg
#perl scripts for counting the number of dirs/files and list in a file.
$DIR="/home/weizg/Documents/data--9.22/rfdata/filter03";
my($num)=@ARGV;
chdir($DIR);
@DIRS = `ls`;
$a=@DIRS-1;
open(NUM,">filels");
print NUM "The num of dirs is $a\n @DIRS:\n";
$total=0;
for($i = 0;$i < @DIRS; $i++)
{
chomp($dir = $DIRS[$i]);
if($num ne ""){$dir=$num;}
if(-d $dir){
chdir($dir);
print "changing to directory $dir\n";
@num2=`ls`;
$b=@num2-2;
print "the num of good files of $dir is $b\n @num2\n";
open(NU,">filels");
print NU "$DIRS[i]";
print NU "the num of good files in $dir is $b\n @num2\n";
print NUM "the num of good files in $dir is $b\n @num2\n";
close(NU);
$totalgood=$b+$totalgood;
chdir(bad1);
@num2b=`ls`;
$c=@num2b-1;
open(NU,">filels");
print NU "$DIRS[i]";
print NU "the num of bad files in $dir is $c\n @num2b\n";
print NUM "the num of bad files in $dir is $c\n @num2b\n";
close(NU);
$totalbad=$c+$totalbad;
chdir("..");
if($dir=$num) {die("The directory $num has finished!\nThe good num is $b;\nThe bad num is $c .\n");}
# IF \n" become  \n " it will output the line of the number.you can try!
}
chdir "..";
}
$total=$totalgood+$totalbad;
print "the num of good total of all is $totalgood\n";
print "the num of bad total of all is $totalbad\n";
print "the num of total of all is $total\n";
print NUM  "the num of good total of all is $totalgood\n";
print NUM  "the num of bad total of all is $totalbad\n";
print NUM  "the num of total of all is $total\n";
close(NUM);
