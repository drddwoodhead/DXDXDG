#!/usr/bin/perl
#initial file to run that runs a spasm based search for matching motifs.

#use strict;
#use warnings;
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/";
use lib "/home/dwoodhead/CaBindData/Modules/";
use runspasmsubsintergrated;
use output;
use input;
use bindingsubs;
use scopsubs;
use clanssubs;
use blastsubs;
use families;



my $familyfile="/home/dwoodhead/CaBind/Input/BCoahuilensispositive.smf";
#my $familyfile="/home/dwoodhead/CaBind/Input/.smf";
#my $familyfile="/home/dwoodhead/CaBind/Input/.smf";
#my $familyfile="/home/dwoodhead/CaBind/Input/.smf";

$familyfile=~m/(^\/.*\/)(.*?)\..*/;
#print "1:$1\n2:$2\n3:$3\n";
$folder=$2;
mkdir "/home/dwoodhead/CaBind/Data/$folder";
$blastdb="/home/dwoodhead/blastdb/nr";#BLAST DATABASE LOCATION
$pdbdir="/home/dwoodhead/CaBindSearch/Data/PDB Database/";
$blastdir="/home/dwoodhead/CaBind/Data/$folder/Blast";
$seqdir="/home/dwoodhead/CaBind/Data/$folder/Sequences/";
$sabledir="/home/dwoodhead/CaBind/Data/$folder/Sable/";
$motiffile=$seqdir."Motiftest.txt";
mkdir $blastdir;
mkdir $seqdir;
mkdir $sabledir;

getfamilydata2($familyfile,$new,$pdbdir,$blastdir,$seqdir,$blastdb,$sabledir);