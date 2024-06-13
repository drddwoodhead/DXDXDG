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


my $file="/home/dwoodhead/CaBindSearch/Resuts/0.5-i-m-b.smh";
my $dbscript="/home/dwoodhead/CaBindSearch/Scripts/mkdb";

my @list=parsexml($file);
my $listcnt=@list;
my $pdblist;
for (my $i=0; $i<$listcnt; $i++)
{
my %hit=%{$list[$i]};
my $id=$hit{pdbid};
$outlist=$outlist."\n/home/dwoodhead/wrkdir/ftpdir/$id.pdb\n";
}
$outlist=$outlist."\n";
`$dbscript $outlist`


