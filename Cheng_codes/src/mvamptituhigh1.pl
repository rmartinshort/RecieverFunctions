############# perl script for moving the files with the amptitude>1 before selecting eqrs
$DIR="/home/weizg/data/NE/f5";    
$dirnew="/home/weizg/data/NE/bad";                                            
chdir($DIR);
@dir=`ls`;
for($j= 0; $j< @dir; $j++)
{
chomp($dir[$j]); 
chdir ($dir[$j]); 
@file=`ls *eqr`;chomp(@file);
for($n=0;$n<@file;$n++)
{
chomp($file[$n]);
print "$file[$n]\n";
$max=`saclst DEPMAX f $file[$n]`;
($temp1,$max) = split(' ',$max); chomp($max);
$min=`saclst DEPMIN f $file[$n]`;
($temp1,$min) = split(' ',$min); chomp($min);
print "$max; $min\n";
if($max>1)
 {
`mv $file[$n] $dirnew `;
 }
}
chdir "..";
}