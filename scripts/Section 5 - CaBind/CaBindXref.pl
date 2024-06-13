#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBind/Modules/";
use lib "/home/dwoodhead/CaBindSearch/Modules";
use input;
use cabindsubs;


#$filelist="/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisASTvDT.results;/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisASTvSVM.results;/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisATGDT.results;/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisATGSVM.results";
#my $rpsfile1="/home/dwoodhead/CaBind/Data/RPS/Bacillus_coahuilensis_scan_Pfam.rpsout";
#my $outfile="/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/TESTBacillus_coahuilensis_Xref.results";

$filelist="/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ASTvDT.results;/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ASTvSVM.results; /home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ATGDT.results;/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ATGSVM.results";
my $rpsfile1="/home/dwoodhead/CaBind/Data/RPS/E_Coil_scan_Pfam.rpsout";
my $outfile="/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110_Xref.results";

my @typeid=('ASTvDT','ASTvSVM','ATGDT','ATGSVM');

my @filelist=split(/;/,$filelist);
my $filecnt=@filelist;
my @currentlist=('null');

for (my $i=0;$i<$filecnt;$i++)
	{
	my $file=$filelist[$i];
	@currentlist=parsecabindresults($file,\@currentlist);
	}
$file=$filelist[0];
@currentlist=parsecabindresults($file,\@currentlist,"N");

@currentlist=rpsblastdata(\@currentlist,$rpsfile1);

print "analysing training data...\n";

$inputfile="/home/dwoodhead/wrkdir/results/families/familygroups-all2.smf";
$rpsfile2="/home/dwoodhead/CaBind/Data/RPS/allfamilygroups.rpsout";

my @familydata=parsexml($inputfile);
@familydata=rpsblastdata2(\@familydata,$rpsfile2);

#@familydata=rpsposdata(\@familydata,$rpsfile);

my %xrefdata=sortrpsdata(\@familydata);

printXrefdata(\@currentlist,$outfile,\@typeid,\%xrefdata);
