#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBind/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules";
use input;
use cabindsubs;


$inputfile="/home/dwoodhead/wrkdir/results/families/familygroups-all2.smf";
$rpsfile2="/home/dwoodhead/CaBind/Data/RPS/allfamilygroups.rpsout";
$outfile="/home/dwoodhead/CaBind/Data/RPS/familyrps.results";

my @familydata=parsexml($inputfile);
@familydata=rpsblastdata2(\@familydata,$rpsfile2);
my %xrefdata=sortrpsdata(\@familydata);

open (OUTFILE,">$outfile");
my @refdataarray=%xrefdata;
my $refdatacnt=@refdataarray;
print "refdataarray:@refdataarray\n";
for (my $i=0;$i<$refdatacnt;$i=$i+2)
{
my $id=$refdataarray[$i];
print OUTFILE "$id\n";
my @rpsref=@{$xrefdata{$id}};
my $rpsrefcnt=@rpsref;
		#print "RPSREF: @rpsref";
		for(my$j=0;$j<$rpsrefcnt;$j++)
			{
			my %hit=%{$rpsref[$j]};
			my $id=$hit{pdbid};
			my $scop=$hit{scop};
			print OUTFILE "			$id	$scop\n";
			}
}