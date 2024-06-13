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

$inputfile="/home/dwoodhead/wrkdir/results/families/familygroups-all2.smf";
$rpsfile="/home/dwoodhead/CaBind/DATA/allfamilygroups.rpsout";

my @familydata=parsexml($inputfile);
@familydata=rpsblastdata2(\@familydata,$rpsfile);
@familydata=rpsposdata(\@familydata,$rpsfile);