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


@hitlist=parsepsout($psscanfile,1);
makegenseqfiles(\@hitlist,$seqdir,$seqfile);