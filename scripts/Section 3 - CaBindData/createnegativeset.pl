#!/usr/bin/perl
# runs spasm just once

use warnings;
use lib "/home/dwoodhead/CaBindData/Modules/";
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/"; 
use runspasmsubsintergrated;
use output;
use psscansubs;
use input;
use scopsubs;
use clanssubs;


my $psscanfile="/home/dwoodhead/wrkdir/scan_out669.out.test";
$negativesetsize="40";#specify approximately how many negative results you want


$runident="negativehits~$negativesetsize";
$resultsdir="/home/dwoodhead/CaBindData/Negative Data/";
$scopdb="/home/dwoodhead/scop/dir.cla.scop.txt_MASTER-pre-scop-2008-12-15.out";

$positivehits="/home/dwoodhead/wrkdir/results/complete/000.5/0.5-i-b-m.smh";
$positivefamily="/home/dwoodhead/wrkdir/results/complete/000.5/0.5-i-b.smf";
#$filename2="/home/dwoodhead/wrkdir/results/negative controls/xxxx:x:0-0.smf";

if(!open(POSITIVE ,$positivehits))
{print "Positive hits file not found\n";
close POSITIVE;
exit;
}
if(!open(NEGATIVE ,$psscanfile))
{print "PS Scan file not found\n";
close NEGATIVE;
exit;
}
print "parsing ps-scan file...\n";
@negativehits=parsepsout($psscanfile); # parses the ps scan file returns initial list of negative resuls




print"parsing positive hits file...\n";
my @positivehits=parsexml($positivehits); # parses the known positive results file

if(open(FAMILY ,$positivefamily))
{
print "...adding family data from file\n";
close FAMILY;
@positivehits=addxmlscopdata(\@positivehits,$positivefamily); #adds scop family data to the positive results file
}
else
{print "no family data file:$familyfile\n";}

$nhitcnt=@negativehits;
print "$nhitcnt Hits\n";

@negativehits=sortnegatives1(\@positivehits,\@negativehits);

#my @negativehits=parsexml($filename2);
$nhitcnt=@negativehits;
print "$nhitcnt Hits\n";

print "adding SCOP classifications...\n";
@negativehits=addscopclass(\@negativehits,"nossm",$scopdb);
@negativehits=sortnegatives2(\@positivehits,\@negativehits);

%hit=%{$positivehits[0]};
$scop = $hit{scop};
$id= $hit{pdbid};
print "ID:$id-SCOP:$scop\n"; 

#$filename2=~m/(.*\/)([^\/]*).smf/ ;
#print "match with : $1 \n";
#print "match with : $2 \n";
$nhitcnt=@negativehits;
print "$nhitcnt Negative hits verified\n";

$nhitcnt=@negativehits;
$filter=$nhitcnt/$negativesetsize;
$filter=int($filter);

$resexpect=3;
@negativehits=filtersuperfamily(\@negativehits,$filter);
$nhitcnt=@negativehits;
print "$nhitcnt Negative hits selected\n";
savehits(\@negativehits,"",$resultsdir,$runident);
savefamilies(\@negativehits,$resultsdir,$runident);
