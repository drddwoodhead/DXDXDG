#!/usr/bin/perl

my $file="/home/dwoodhead/CaBindTools/Input/psipredlistmore2.txt";
my $outfile="/home/dwoodhead/CaBindTools/Input/psipredlistmore2.smh";

open(FILE,"<$file");
open(OUTFILE,">$outfile");
print OUTFILE "<XML>\n<smh type='bhits'>\n";
my @array=<FILE>;
for $line(@array)
{
if ($line=~m/(\d*):(-):(\d{1,4})-(\d{1,4})/)
{
my $mid=$3+2;
print OUTFILE "<hit id=$1:$2:$3-$4 run=>\n";
print OUTFILE "<file>$1<\\file>\n";
print OUTFILE "<chain>$2<\\chain>\n";
print OUTFILE "<res>$3<\\res>\n";
print OUTFILE "<res>$mid<\\res>\n";
print OUTFILE "<res>$4<\\res>\n";
print OUTFILE "<\\hit>\n\n";
}
else
{}
}
print OUTFILE "<\smh>\n<\XML>\n";
close FILE;
close OUTFILE;