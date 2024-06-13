#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBindData/Modules/";
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/";
use scopsubs;
use output;
use blastsubs;
use input;
use families;
use dtdata;
use svmdata;
use testcases;
use cabindsubs;

#$seqfile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis_test.txt";
#$seqfile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis.txt";
#$seqfile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110_test.txt";
$seqfile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110.txt";
#$dtrootfilename="/home/dwoodhead/CaBind/Input/ASTv/ASTv";
#$dtrootfilename="/home/dwoodhead/CaBind/Input/HYC4/HYC4";
$dtrootfilename="/home/dwoodhead/CaBind/Input/ATG/ATG";
#$svmmodelfile="/home/dwoodhead/CaBind/Input/ASTv/ASTv.model";
#$svmmodelfile="/home/dwoodhead/CaBind/Input/HYC4/HYC4.model";
$svmmodelfile="/home/dwoodhead/CaBind/Input/ATG/ATG.model";


#$genomefile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis.ref";
#$rpsfile="/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensis_scan_Pfam.rpsout";

$genomefile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110.ref";
$rpsfile="/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110_scan_Pfam.rpsout";

$dtrootfilename=~m/^.*\/(.+)/;
my $exp=$1;
#my $exp="HYT";
#my $exp="ATG";

my $expdt=$exp."DT";
my $expsvm=$exp."SVM";

$seqfile=~m/(^\/.*\/)(.*?)\..*/;
#print "1:$1\n2:$2\n3:$3\n";
$runident=$2;

my $iddata=0;
my $blast; #used to ensure the data is only collected if needed.
my $aasize="0"; #AT Threshold size L/S as percentage across blast, 0;off, 1-14;threshold value, v;different threshold for each residue.
my @aasizearray=(4,3,11,4,4,8,3,8,9,3,2,9,1);#when $aasize="v" use to specify threshold for each residue position.
my $aatype=1;
my $aahyd="0";	#HYC Threshold hydrophobicity Chacracter HYPHO/HYPIL as percentage across blast.
my @aahydarray=(0,0,0,0,0,0,0,0,0,0,0,0,0);
my $ss=0;	#SSX Secondary structure data
my $con=0;	#CON Conserved residues
my $sol=0;	#SOL Solvent accesability
my $blast=1;
my $surround=4;

$blastdir="/home/dwoodhead/wrkdir/blast/";
$seqdir="/home/dwoodhead/wrkdir/sequences/";
$motiffile=$seqdir."Motiftest.txt";
$blastdb="/home/dwoodhead/blastdb/nr";
$resultsdir="/home/dwoodhead/CaBind/Data/$runident/";
$blastdir="/home/dwoodhead/CaBind/Data/$runident/Blast/";
$seqdir="/home/dwoodhead/CaBind/Data/$runident/Sequences/";
$sabledir="/home/dwoodhead/CaBind/Data/$runident/Sable/";
$psipreddir="/home/dwoodhead/CaBind/Data/$runident/Psipred/";

$hitfile=$resultsdir.$runident.".smh";

my @hitlist=parsexml($hitfile);

my @families=sortbyclass(\@hitlist);

savefamilies(\@families,$resultsdir,$runident);
$familyfile=$resultsdir.$runident.".smf";

@familydata=processfamilydata($familyfile,$surround,$aasize,\@aasizearray,$aahyd,\@aahydarray,$con,$sol,$ss,$blast,$blastdir,$psipreddir,$sabledir);

@familydata=genomedata(\@familydata,$genomefile);
@familydata=rpsblastdata(\@familydata,$rpsfile);

%dtresults=classifybydt($dtrootfilename,\@familydata,$aasize,\@aasizearray);
saveclassresults(\%dtresults,$familyfile,$dataid,$resultsdir,$expdt);
savegenresults(\%dtresults,$familyfile,$dataid,$resultsdir,$expdt);

%svmresults=classifybysvm($svmmodelfile,\@familydata,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runident,$aa,$size,$abssdata,$hyd,$res,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,\@aatypearray,$resultsdir,$runid);
saveclassresults(\%svmresults,$familyfile,$dataid,$resultsdir,$expsvm);
savegenresults(\%svmresults,$familyfile,$dataid,$resultsdir,$expsvm);
