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


my $familyfile="/home/dwoodhead/CaBindTrain/Input/testdata.smf";

$blastdb="/home/dwoodhead/blastdb/nr";#BLAST DATABASE LOCATION
$pdbdir="/home/dwoodhead/CaBindSearch/Data/PDB Database/";
$blastdir="/home/dwoodhead/CaBindData/Data/Blast/";
$seqdir="/home/dwoodhead/CaBindData/Data/Sequences/";
$sabledir="/home/dwoodhead/CaBindData/Data/Sable/";
$motiffile=$seqdir."Motiftest.txt";


getfamilydata($familyfile,,$pdbdir,$blastdir,$seqdir,$blastdb,$sabledir);