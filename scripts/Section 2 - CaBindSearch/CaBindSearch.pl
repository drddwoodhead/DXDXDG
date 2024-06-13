#!/usr/bin/perl
#initial file to run that runs a spasm based search for matching motifs.

#use strict;
#use warnings;
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
$rmscutoff="1";
$acmiss="2";
$scmiss="2.5";
#$acmiss=$rmscuttoff+1;
#$scmiss=$rmscuttoff+1;
$spasmDB="/home/dwoodhead/spasm/db/spasm_june_2010.lib";
$scopdb="/home/dwoodhead/scop/dir.cla.scop.txt_MASTER-pre-scop-2008-12-15.out";

$startfile="/home/dwoodhead/CaBindSearch/Input/1exr_A_20-24.pdb";


if(!open(STARTFILE ,$startfile))
{
print "Starting PDB file not found\n $startfile\n";
close STARTFILE;
exit;
}

#@results=`$spasm $database $pdb $rmscutoff $inid \"$substitutions\"

#########AUTO INPUTS#######
$forcespasm='0';
$resnumber=3;
$resexpect=3;

$spasmDB=~m/(^.*\/)((.+)\.lib)/;
$sdb=$3;
$spasmdatabase=$2;
$spasmdir=$1;
#print "1:$1,2:$2,3:$3\n";
$scopdir='/home/dwoodhead/scop/';
$workdir='/home/dwoodhead/CaBindSearch/'; 
$macro=$workdir.'temp.omac';

$startfile=~m/(^.*\/)((.+)\.pdb)/;
$seedpdb=$3;
$seedfile=$2;
$seeddir=$1;
#print "1:$1,2:$2,3:$3\n";
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

$runident=$seedpdb;

#$ftpscript=$workdir.'Script/spasmftp';# location of script that downloads pdb files.
$spasm=$workdir.'Scripts/DXDXDGg'; #### location of shell script to run spasm

########################################################################


$absolutehits=0;
@runlist;
@masterhitarray;
$seedpdb=~s/_/:/g;
@mastertrylist=($seedpdb);




                
@masterhitarray=initialisearray($startfile,$sitedir);
my %each=%{$masterhitarray[0]};
my $ref=\@{$each{resnum}} ;
#print "making pdb ".$each{pdbid}." residues "."@{$ref}"."\n";
mkpdb ($each{pdbid}, $ftpdir, $ref, $each{chain},$sitedir);
@masterhitarray=addmetalbinding(\@masterhitarray,$ftpdir,$sitedir);
savehits(\@masterhitarray,"",$resultsdir,$runident);


$runlistcount=@runlist;
$trylistcount=@mastertrylist;
print "runlist count = $runlistcount, trylist count = $trylistcount\n\n";
my$iteration=0;
while ($trylistcount>$runlistcount)
#for(my $x=0;$x<1;$x++)
{
$totalruns=@runlist;
#print "$totalruns runs so far\n";
#print "sitedir:$sitedir\n";
print "***********************Iteration $iteration***********************\n";
%lists=iterativerun(\@mastertrylist, $spasmDB, $ftpdir, $sitedir,$rmscutoff,$acmiss,$scmiss,$logdir,$resultsdir,$seedpdb,$runlistcount,\@runlist,$iteration);
@mastertrylist=@{$lists{trylist}};
@runlist=@{$lists{runlist}};
$runlistcount=@runlist;
$trylistcount=@mastertrylist;
print "Hits tried $runlistcount, Hits in list $trylistcount\n\n";
$iteration++;
}
 
my $filename=$resultsdir.$runident."\.smh";
 # open(OLDFILE,"<$filename");
  @oldresults=<OLDFILE>;
  $line1=pop @oldresults;
  $line2=pop @oldresults;
  $line3=pop @oldresults;
  $line4=pop @oldresults;  
#  open(RESULTS,">$filename");
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

 print    "********************************************************************************\n SPASM Analysis complete\n********************************************************************************";

@hitlist=parsexml($filename);


@hitlist=addscopclass(\@hitlist,"nossm",$scopdb);

my @families=sortbyclass(\@hitlist);

savefamilies(\@families,$resultsdir,$runident);

