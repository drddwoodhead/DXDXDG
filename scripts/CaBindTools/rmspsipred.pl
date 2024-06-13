#!/usr/bin/perl
#RMS
use warnings;
use lib "/home/dwoodhead/CaBindTools/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBindData/Modules/";
use input;
use blastsubs;
use rmssubs;
#use clanssubs;
use bindingsubs;
use output;


$id1="1exr:A:24-29";
#$rmsfilename="/home/dwoodhead/CaBindTools/Input/rmslist.smh";
$predfilename="/home/dwoodhead/CaBindTools/Input/psipredlistmore2.smh";
#$rmsoutfilename="/home/dwoodhead/CaBindTools/Input/rms.out";
#$rmsatomtypes="TOOL";

$workdir='/home/dwoodhead/CaBindTools/'; 
$ftpdir=$workdir.'Input/db/';
#$sitedir=$workdir.'Input/sitepdbs/';
$outdir=$workdir.'Input/db/';
$blastdir=$workdir.'Input/blastdir/';
$rmsatomtypes='TOOL';
$id1='1EXR:A:20-24';

$predfilename=~m/(.*\/)([^\/]*).smh/ ;
print "match with : $1 \n";
print "match with : $2 \n";
 
$resultsdir=$1;
$runident=$2;
#$resexpect=3;

#print"parsing rmslist file...\n";
#my@array=parsexml($rmsfilename);
#my $hitlistcount=@array;
#open(RMSOUTFILE,">$rmsoutfilename");
#print "hitlistcount=$hitlistcount";
#for (my $i=0; $i<$hitlistcount; $i++)
	#{
         #my %hit=%{$array[$i]};
	 #my $id=$hit{file};
	 #my $pdbid=$hit{pdbid};
	 #my $seqfile="$ftpdir/$pdbid.fas";
	 #my $chain=$hit{chain};
	 # #print "CHAIN $chain\n";
	 #my @res=@{$hit{resnum}};
	 # print "mkpdb($pdbid,$ftpdir,\@res,$chain,$sitedir,1);\n";
	 #mkpdb($pdbid,$ftpdir,\@res,$chain,$sitedir,"1");
	 #print "runlsqman($id1,$id,$rmsatomtypes);\n";
	 #my $rms=runlsqman($id1,$id,$rmsatomtypes,$sitedir);
	 #print RMSOUTFILE "$id,$id1,$rms\n";
	#}
print"parsing psipredlist file...\n";
my@array2=parsexml($predfilename);
my $hitlistcount=@array2;


for (my $i=0; $i<$hitlistcount; $i++)
	{
        my %hit=%{$array2[$i]};
	 my $id=$hit{file};
	 my $pdbid=$hit{pdbid};
	 my $seqfile="$ftpdir$pdbid.fas";
	 print "seqfile:$seqfile\n";
	 my $chain=$hit{chain};
	 my @res=@{$hit{resnum}};
	  
	 psipred(\%hit, $range, $blastdir, $outdir,$seqfile);
	}
