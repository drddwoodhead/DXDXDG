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
$outputfile="/home/dwoodhead/CaBind/Data/allfamilygroups.fas";

$seqdir="/home/dwoodhead/wrkdir/sequences/";
$ftpdir="/home/dwoodhead/wrkdir/ftpdir/";

open(OUTFILE,">$outputfile");

my @hitlist=parsexml($inputfile);

my $hitcnt=@hitlist;

for(my$i=0;$i<$hitcnt;$i++)
{
my %hit=%{$hitlist[$i]};
my $scop=$hit{scop};
my @pdblist=@{$hit{pdblist}};
my $pdbcnt=@pdblist;
#print @pdblist;

for (my$j=0;$j<$pdbcnt;$j++)
{

my$id=$pdblist[$j];
#print "ID:$id\n";
print OUTFILE ">$id\n";

$id=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
#print "$1,$2,$3,$4,$5,$6,$7,$8,$9\n";
my $pdbid=$1;
my $chain=$2;
my @res=($3,$3+2,$4);
#print "ID:$pdbid,CHAIN:$chain,RES1:$res[0],RES3:$res[2]\n";

my %newhit=('resnum',\@res,'chain',$chain,'pdbid',$pdbid);
my $seq=getseq(\%newhit,"","",$ftpdir);
print OUTFILE "$seq"; 
}
} 