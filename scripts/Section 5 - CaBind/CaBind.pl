#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBindData/Modules/";
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/";

use output;
use blastsubs;
use input;
use families;
use runspasmsubsintergrated;
use output;
use psscansubs;

my $phiblast=1;
my $psipred=0;
my $sable=0;



$runident="xxxx:x:0-0";
#my $seqfile="/home/dwoodhead/CaBind/Input/test.txt";
my $seqfile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110.txt";
#my $seqfile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110_test.txt";
#my $seqfile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis_test.txt";
#my $seqfile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis.txt";


$seqfile=~m/(^\/.*\/)(.*?)\..*/;
#print "1:$1\n2:$2\n3:$3\n";
$runident=$2;
$resultsdir="/home/dwoodhead/CaBind/Data/$runident/";
`mkdir \"$resultsdir\"`;
my $psscanfile=$resultsdir.$runident.".scn";
#my $familyfile=$resultsdir.$runident.".smh";
#print "RUNID:$runident\n";
#print "RESULTDIR:$resultsdir\n";
#print "PSSCANFILE:$psscanfile\n";
#print "FAMILYFILE:$familyfile\n";

$blastdir="/home/dwoodhead/CaBind/Data/$runident/Blast/";
$seqdir="/home/dwoodhead/CaBind/Data/$runident/Sequences/";
$sabledir="/home/dwoodhead/CaBind/Data/$runident/Sable/";
$psipreddir="/home/dwoodhead/CaBind/Data/$runident/Psipred/";

$motiffile=$seqdir."Motiftest.txt";
`mkdir \"$blastdir\"`;
`mkdir \"$seqdir\"`;
`mkdir \"$sabledir\"`;
`mkdir \"$psipreddir\"`;
$blastdb="/home/dwoodhead/blastdb/nr";
chdir($psipreddir);

print "Running ps scan to search for $pattern motif...\n";
$pattern="D-x-[DNS]-x[DNS]";
runpsscan($seqfile,$psscanfile,$pattern); # runs ps scan and on $seqfile saves results to $psscanfile
@hitlist=parsepsout($psscanfile,1); #parses $psscanfile and puts results into an array:- @hitlist;
savehits(\@hitlist,"",$resultsdir,$runident); #creates an XML file from the array @hitlist

makegenseqfiles(\@hitlist,$seqdir,$seqfile);

#print "Retiveing data associated with results of ps scan...\n";

getgendata(\@hitlist,"",$seqdir,$blastdir,$sabledir,$phiblast,$psipred,$sable); 
saveclassresults(\%dtresults,$familyfile,$dataid,$resdir);
