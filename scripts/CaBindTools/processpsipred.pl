#!/usr/bin/perl

my $file="/home/dwoodhead/CaBindTools/Input/psipredlistmore.txt";


open(FILE,"<$file");
my @array=<FILE>;
my $outfileall="/home/dwoodhead/CaBindTools/Input/blastdir/all.jaa";
open(OUTFILEALL,">$outfileall");
print OUTFILEALL "JALVIEW_ANNOTATION\n\n";

for $line(@array)
{
$line=~m/(.*?):.*/ ;
print "ID:$1\n";
my $id=$1;
my $infile="/home/dwoodhead/CaBindTools/Input/blastdir/$id.psipred";
my $outfile="/home/dwoodhead/CaBindTools/Input/blastdir/$id.jaa";
open (INFILE,"<$infile");
my @inarray=<INFILE>;
open(OUTFILE,">$outfile");
print OUTFILE "JALVIEW_ANNOTATION\n\n";
my $sstring="NO_GRAPH	$id	Secondary Structure	";
my $cstring="NO_GRAPH	conf $id	Confidence	";
for $inline(@inarray)
{
if ($inline=~m/^Pred\:\s(.*)\s/)
{
#print "LINE:  $1\n";
my $string=join('|', split(//, $1));
#print "B4:$string\n";
$string=~s/C//g;
#print "AF:$string\n";
$sstring=$sstring."$string|";
}
elsif ($inline=~m/^Conf\:\s(.*)\s/)
{
#print "LINE:  $1\n";
my $string=join('|', split(//, $1));
#print "CString:$string\n";
$cstring=$cstring."$string|";
}
else {}

}
print OUTFILE $sstring."\nCOLOUR	        $id	000000\n";
print OUTFILEALL $sstring."\nCOLOUR	$id	000000\n\n";
print OUTFILE $cstring."\nCOLOUR	        conf $id	000000\n";
print OUTFILEALL $cstring."\nCOLOUR	conf $id	000000\n\n";
close INFILE;
close OUTFILE;
}

close OUTFILEALL;
close FILE;
