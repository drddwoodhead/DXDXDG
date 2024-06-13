#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBindData/Modules/";
use lib "/home/dwoodhead/CaBindTrain/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules/";
use lib "/home/dwoodhead/CaBind/Modules/";
use runspasmsubsintergrated;
use output;
use blastsubs;
use input;
use families;
use dtdata;
use testcases;

my $iddata=0;

my $aasize="v"; #AST Threshold size L/S as percentage across blast, 0;off, 1-14;threshold value, v;different threshold for each residue.
my @aasizearray=(4,3,11,4,4,8,3,8,9,3,2,9,1);#when $aasize="v" use to specify threshold for each residue position.
my $aatype=0;	#ATG Amino acid group as percentage +,-,P,N,G
my $aahyd="0";	#HYC Threshold hydrophobicity Chacracter HYPHO/HYPIL as percentage across blast.
my @aahydarray=(0,0,0,0,0,0,0,0,0,0,0,0,0);

my $ss=0;	#SSX Secondary structure data
my $con=0;	#CON Conserved residues
my $sol=0;	#SOL Solvent accesability


my $abssdata=0; #ABS Threshold size of hit only L/S, 0;off, 1;on
my $size="0";	#ABW Amino acid weight of hit only 
my @sizearray=(0,0,0,0,0,0,0,0,0,0,0,0,0);
my $abstdata=0; #ABG Amino acid group of hit only +,-,P,N,G 0;off, 1;on
my $aa=0;	#ABR Amino Acid one letter code
my $hyd=0;	#ABH Amino acid hydrophobicity of hit only
my $abshdata=0; #ABC threshold hydrophobicity character of hit only HYPHO/HYDPHIL, 0;off, 1;on

my $avesize=0;  #AWA Average amino acid Weight
my @avesizearray=(0,0,0,0,0,0,0,0,0,0,0,0,0);
my $aaper=0;	#ATR Amino acid residue as percentage
my $avehyd=0;	#HYA Average hydrophobicity					
my @avehydarray=(0,0,0,0,0,0,0,0,0,0,0,0,0);

#my $experiment="test";
#my $experiment="all attrib";
my $experiment="selective attrib";
#my $experiment="variable threshold";
#my $experiment="AA size threshold or type";
#my $experiment="AAsize Threshold";
#my $experiment="AAhyd Threshold";
#my $experiment="AA combined";
#my $experiment="threshold and continuous";
#my $experiment="best threshold";
my $surround=4;

@trackskip;

my $familyfile="/home/dwoodhead/CaBindTrain/Input/testdata.smf";


$familyfile=~m/^.*\/(.+)\.smf/;
my $fileid=$1;

print "match with : $1 \n";
print "match with : $2 \n";
$blastdir="/home/dwoodhead/wrkdir/blast/";
$seqdir="/home/dwoodhead/wrkdir/sequences/";
$motiffile=$seqdir."Motiftest.txt";

$blastdb="/home/dwoodhead/blastdb/nr";
my $blast; #used to ensure the data is only collected if needed.
if($abstdata==0 && $abssdata==0 && $aasizedata==0 && $aatype==0 && $aa==0 && $size==0 && $hyd==0 && $size ne "v" &&$size !=0 && $aasize ne "v" &&  $aasize!=0 && $aahyd==0 && $avesize==0 && $avehyd==0)
{
#print "no blast\n";
$blast=0;
}
else
{
$blast=1;
}
#print "BLAST==$blast\n"; 
my $raasize; #variable used for folder/file naming when $aasize="v"
for (my $i=0;$i<14;$i++)
	{
	my $each=$aasizearray[$i];
	if($each!=0)
		{
		my $resnum=$i-4;
		$raasize=$raasize."$each"; #filenaming R5:6 means Residue position 5 has a threshold value of 6
		}	
	}
my $rsize; 
for (my $i=0;$i<14;$i++)
	{
	my $each=$sizearray[$i];
	if($each!=0)
		{
		my $resnum=$i-4;
		$rsize=$rsize."$resnum,";
		}	
	}
my $raahyd; #variable used for folder/file naming when $aasize="v"
for (my $i=0;$i<14;$i++)
	{
	my $each=$aahydarray[$i];
	if($each!=0)
		{
		my $resnum=$i-4;
		$raahyd=$raahyd."$each"; #filenaming R5:6 means Residue position 5 has a threshold value of 6
		}	
	}
print "processing family data...\n";
@familydata=processfamilydata($familyfile,$surround,$aasize,\@aasizearray,$aahyd,\@aahydarray,$con,$sol,$ss,$blast,$blastdir);
$famcnt=@familydata;
print "FAMILYDATA:$famcnt\n";
@families=pickpositive($familyfile);
$pzcnt=@families;
print "posative Families:$pzcnt\n";

#### this Section is concerened with giving unique names to differnt sets of variables.
my $runident=$fileid."_";
if($abstdata==1){$runident=$runident."ABS_";}
if($abssdata!=0){$runident=$runident."ABT_";}
if($aasize!=0){$runident=$runident."ASTR$aasize"."_";}
#if($aasize eq "v"){$runident=$runident."AST$raasize"."_";}
if($aasize eq "v"){$runident=$runident."ASTv"."_";}
if($aatype==1){$runident=$runident."ATG_";}
if($ss==1){$runident=$runident."SSX_";}
if($con==1){$runident=$runident."CON_";}
if($sol==1){$runident=$runident."SOL_";}
if($aa==1){$runident=$runident."ABR_";}
if($size eq "v"){$runident=$runident."ABWR$rsize"."_";}
if($size==1  && $res==0){$runident=$runident."ABW_";}
if($size==1 && $res!=0){my $myres=$res-5;
			$runident=$runident."ABWR$myres"."_";}
if($hyd!=0){$runident=$runident."ABH_";}
if($aahyd!=0){$runident=$runident."HYC$aahyd"."_";}
#if($aahyd eq "v"){$runident=$runident."HYC$raahyd"."_";}
if($aahyd eq "v"){$runident=$runident."HYCv"."_";}
if($avesize==1){$runident=$runident."AWA_";}
if($avehyd==1){$runident=$runident."HYA_";}


$dir="/home/dwoodhead/CaBindTrain/Data/";
$dtresultsdir=$dir."DT_".$runident."/";
`mkdir $dir`;
`mkdir "$dtresultsdir"`;
my $family=@families;
$family=$family+1;

my $errorfile="$dtresultsdir$runident.errors";
open (ERRORFILE, ">$errorfile");
print ERRORFILE "$runident\n";
print ERRORFILE "Skip\t \% Errors\tMissclassified\tErrorlist\n";
close ERRORFILE;

my $statsfile="$dtresultsdir$runident.stats";
open (STATSFILE, ">$statsfile");
print STATSFILE "$runident\n";
print STATSFILE "Skip\tTrue Positive\tTrue Negative\tFalse Positive\tFalse Negative\tErrors\n";
close STATSFILE;


for (my $i=0;$i<$family;$i++)
{
my $skip=$i;
#my $runid=$runident."skip$skip";
my $runid="skip$skip";

my $resdir=$dtresultsdir."skip$skip/";
`mkdir "$resdir"`;

my %skipfamily=%{$families[$skip]};
$skipfamily=$skipfamily{scop};
print "Skip:$skip:$skipfamily\n";
if ($skipfamily eq "")
{
$skipfamily="none";
}

my @familydataskip=skip(\@familydata,$skipfamily);
saveDTdata(\@familydataskip,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runid,$aa,$size,$abssdata,$hyd,$res,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,$avesize,$avehyd,$resdir,$runid);
saveDTnames($aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runid,$aa,$size,$abssdata,$hyd,$res,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,$avesize,$avehyd,$resdir,$runid);
my $track=@familydataskip;
my $tot=@familydata;
#print "\nSkip $skipfamily $track/$tot cases used\n";
push (@trackskip,$track);

$dtrootfilename=$resdir.$runid;
$rootident=$resdir.$runid;
chdir $dtresultsdir;
print "/home/dwoodhead/C4.5/R8/Src/c4.5 -f \"$rootident\">\"$dtrootfilename.dt\"\n";
$any1=`/home/dwoodhead/C4.5/R8/Src/c4.5 -f \"$rootident\">\"$dtrootfilename.dt\"`;
print "/home/dwoodhead/C4.5/R8/Src/c4.5 -f \"$rootident\">\"$dtrootfilename.r\"\n";
$any2=`/home/dwoodhead/C4.5/R8/Src/c4.5rules -f \"$rootident\">\"$dtrootfilename.r\"`;

print "Classify by DTs\n\n";

%dtresults=classifybydt($dtrootfilename,\@familydata,$aasize,\@aasizearray);
%dtresults2=classifybydt($dtrootfilename,\@familydataskip,$aasize,\@aasizearray);
saveclassresults(\%dtresults,$familyfile,$dataid,$resdir);
saveclassresults(\%dtresults2,$familyfile,$dataid,$resdir);
%errors=getDTerrors("$dtrootfilename.dt",\@familydata,\%dtresults,\%dtresults2,$runident,$dtresultsdir,$skip);
}
my $errorfile="$dtresultsdir$runident.errors";
my $treefile="$dtresultsdir$runident.trees";
my $statsfile="$dtresultsdir$runident.stats";
my $mergefile="$dtresultsdir$runident.merge";
mergefile($errorfile,$treefile,$statsfile,$mergefile);
print "$errorfile,$statsfile\n";
#print "@trackskip";



