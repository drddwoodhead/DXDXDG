#!/usr/bin/perl
# runs spasm just once

#use warnings;
use lib "/home/dwoodhead/CaBind/Modules/";

use cabindsubs;


my $filelist="/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisASTvDT.results;/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisASTvSVM.results; /home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisATGDT.results;/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensisATGSVM.results";
$genomefile="/home/dwoodhead/CaBind/Input/Bacillus_coahuilensis.ref";

#$filelist="/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ASTvDT.results;/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ASTvSVM.results; /home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ATGDT.results;/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110ATGSVM.results";
#$genomefile="/home/dwoodhead/CaBind/Input/E_Coli_K12_W3110.ref";

my $outfile="/home/dwoodhead/CaBind/Data/Bacillus_coahuilensis/Bacillus_coahuilensis_ALL2.results";
#my $outfile="/home/dwoodhead/CaBind/Data/E_Coli_K12_W3110/E_Coli_K12_W3110_ALL2.results";
#my $outfile="/media/Data/Documents/PhD/Results Section 5/ALL.results";

#my @typeid=('ASTvDT','ASTvSVM','ATGDT','ATGSVM','HYC4DT','HYC4SVM');
my @typeid=('ASTvDT','ASTvSVM','ATGDT','ATGSVM',);

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
@currentlist=genomedata(\@currentlist,$genomefile);

printcollateddata(\@currentlist,$outfile,\@typeid);
