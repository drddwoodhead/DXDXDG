#!/usr/bin/perl
#initial file to run that runs a spasm based search for matching motifs.

#use strict;
use warnings;
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/";
use runspasmsubsintergrated;
use output;
use input;
use bindingsubs;
use scopsubs;
use clanssubs;


#######SPASM INPUTS#######
$rmscutoff=1;
$acmiss="\n";
$scmiss="\n";
#$acmiss=$rmscuttoff+1;
#$scmiss=$rmscuttoff+1;
$spasmDB="/home/dwoodhead/spasm/spasm_100_jul08.lib";

$startfile="/home/dwoodhead/CaBindTrain/Input/1exr_A_20-24.pdb";
$resumefilename="";


if(!open(STARTFILE ,$startfile))
{
print "Starting PDB file not found\n $startfile\n";
close STARTFILE;
exit;
}

#@results=`$spasm $database $pdb $rmscutoff $inid \"$substitutions\"

#########AUTO INPUTS#######
$iterations='i'; # no longer specifys number of iterations but [i]nfinate, for file naming
$forcespasm='0';
$resnumber=3;
$resexpect=3;

$spasmDB=~m/(^.*\/)((.+)\.lib)/;
$sdb=$3;
$spasmdatabase=$2;
$spasmdir=$1;
print "1:$1,2:$2,3:$3\n";
$scopdir='/home/dwoodhead/scop/';
$workdir='/home/dwoodhead/CaBindTrain/'; 
$macro=$workdir.'temp.omac';

$startfile=~m/(^.*\/)((.+)\.pdb)/;
$seedpdb=$3;
$seedfile=$2;
$seeddir=$1;
print "1:$1,2:$2,3:$3\n";
my $resultsdir=$workdir.'Results/'.$sdb.'/'.$seedpdb.'/'.$rmscutoff.'/';
my $ftpdir=$workdir.'Data/PDB Database/';
my $sitedir=$workdir.'Data/Motif Files/';
my $logdir=$workdir."Data/Logs/$sdb/$rmscutoff/";
mkdir $workdir.'Data/';
mkdir $workdir.'Results/';
mkdir $workdir.'Results/'.$sdb.'/';
mkdir $workdir.'Results/'.$sdb.'/'.$seedpdb.'/';
mkdir $resultsdir;
mkdir $ftpdir;
mkdir $sitedir;
mkdir $workdir.'Logs/';
mkdir $workdir."Logs/$sdb/";
mkdir $logdir;
mkdir $workdir.'Results/';

$runident=$seedpdb."\-".$rmscutoff."\-".$iterations;

#$ftpscript=$workdir.'Script/spasmftp';# location of script that downloads pdb files.
$spasm=$workdir.'Scripts/DXDXDGg'; #### location of shell script to run spasm

########################################################################


my %retrivedresults=parsexml($resumefilename);

my @mastertrylist=@{$retrivedresults{hitarray}};
my $runcount=$retrivedresults{runcount};


$runlistcount=@runlist-1;
$trylistcount=@mastertrylist;
print "initial runlist count = $runlistcount, trylist count = $trylistcount\n\n";


while ($trylistcount>$runlistcount)
{
$totalruns=@runlist;
print "$totalruns runs so far\n";
@mastertrylist=iterativerun(\@mastertrylist, $spasmDB, $ftpdir, $sitedir,$rmscutoff,$logdir,$resultsdir,$seedpdb);
$runlistcount=@runlist-1;
$trylistcount=@mastertrylist;
print "runlist count = $runlistcount, trylist count = $trylistcount\n\n";
$x++;
}

my $filename=$resultsdir.$runident."\.smh";
  open(OLDFILE,"<$filename");
  @oldresults=<OLDFILE>;
  $line1=pop @oldresults;
  $line2=pop @oldresults;
  $line3=pop @oldresults;
  $line4=pop @oldresults;  
  open(RESULTS,">$filename");
  for $line (@oldresults)
    {
    print RESULTS "$line";
    }
  print RESULTS "$line3";
  print RESULTS "$line2";
  print RESULTS "$line1";
  close RESULTS;
  close OLDFILE;
 
 $finalhits=@mastertrylist;
 print ("\n Results\n");
 print "Total number of hits $absolutehits\n";
 print "Non redundant hits $finalhits\n";

 print    "********************************************************************************\n Analysis complete\n********************************************************************************";




@mastertrylist=addscopclass(\@mastertrylist,"nossm");

my @families=sortbyclass(\@array);

savefamilies(\@families);