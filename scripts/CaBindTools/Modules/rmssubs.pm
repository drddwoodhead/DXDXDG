sub runlsqman
{
$ID1=$_[0];
$ID2=$_[1];
$atomtypes=$_[2];
$pdbdb=$_[3];

checksitepdb($ID1,$pdbdb);
checksitepdb($ID2,$pdbdb);

$ID1=~m/(.*?):(.):(\d{1,4})-(\d{1,4})/;
$chain1="A";
$start1=$3;
$end1=$4;
$fileID1="$1_$2_$3-$4";

$ID2=~m/(.*?):(.):(\d{1,4})-(\d{1,4})/;
$chain2="A";
$start2=$3;
$end2=$4;
$fileID2="$1_$2_$3-$4";


print "range of first:$chain1:$start1-$end1 range of second:$chain2:$start2-$end2 \n";

if ($atomtypes eq "DE")
{
$lsqmanscript="/home/dwoodhead/wrkdir/lsqmanscriptDE";
}
elsif ($atomtypes eq "CA")
{
$lsqmanscript="/home/dwoodhead/wrkdir/lsqmanscriptCA";
}
elsif($atomtypes eq "TOOL")
{
$lsqmanscript="/home/dwoodhead/CaBindTools/Scripts/lsqmanscriptTOOL"
}
my @lsqresults=`$lsqmanscript $fileID1 $fileID2 $chain1 $start1 $end1 $chain2 $start2 $end2 $pdbdb`;
print "RESULT of SCRIPT: @lsqresults\n";
my $rms;

#my $regex=".{10,15}atoms have an RMS distance of\s*\d{1,2}\.\d{3}\s* A";
for $line(@lsqresults)
{
if ($line =~m/^.*atoms have an RMS distance of(.*)A/)
          {
	#print "$line\n";
        $rms=$1;
	$rms=~s/\s//g;
          }

}
return $rms;

}

sub makermshash
{
@hitlist=@{$_[0]};
#my $atoms=$_[1];
my $hitlistcount=@hitlist;

my %rmshash;

my $pdbdb="/home/dwoodhead/wrkdir/sitepdb/";
for (my $i=0; $i<$hitlistcnt; $i++)
	{
	for (my $j=$i; $j<$hitlistcnt; $j++)
	       {
		my %hit1=%{$hitlist[$i]};
		my %hit2=%{$hitlist[$j]};
		
		my $ID1=$hit1{file};
		my $ID2=$hit2{file};
		
		my $rms=runlsqman($ID1,$ID2, $rmsatomtypes,$pdbpd);

		if($rmshash{$ID1})
		{
		#print "$ID1 found in hash\n";
		my %ID1HS=%{$rmshash{$ID1}};
		$ID1HS{$ID2}=$rms;
		$rmshash{$ID1}=\%ID1HS;
		}
		else
		{
		#print "$ID1 not found in hash\n";
		my %ID1HS;
		$ID1HS{$ID2}=$rms;
		$rmshash{$ID1}=\%ID1HS;
		}

		if($rmshash{$ID2})
		{
		#print "$ID2 found in hash\n";
		my %ID2HS=%{$rmshash{$ID2}};
		$ID2HS{$ID1}=$rms;
		$rmshash{$ID2}=\%ID2HS;
		}
		else
		{
		#print "$ID2 not found in hash\n";
		my %ID2HS;
		$ID2HS{$ID1}=$rms;
		$rmshash{$ID2}=\%ID2HS;
		}
	       }
	}
return %rmshash;
}

sub savermsnotworking
{
%rmstable=%{$_[0]};
@rmstab=@{%rmstable};

my $filename=$resultsdir.$runident."\.smr";
open(MATRIXFILE,">$filename"); # open new sequence file


print MATRIXFILE  "<XML>\n<smh type='rms'>\n";

for (my $i=0; $i<$hitlistcount; $i++)
	{
	my$ID1=
	print MATRIXFILE "<RMS ID=".$ID1." ID=".$ID2 ."> ".$rms." <\\RMS>\n";
	}
print MATRIXFILE "<\\smh><\\XML>\n";
close MATRIXFILE;
}

sub addrmsdata
{
 my @hitlist=@{$_[0]};
 my $filename=${$_[1]};

 my %rmslist=parsexml(\$filename);

my $hitlistcnt=@hitlist;

for (my $i=0; $i<$hitlistcnt; $i++)

{
my %hit=%{$hitlist[$i]};
my $hitid=$hit{file};
#print"hit = $hitid\n";
my%rmsdata=%{$rmslist{$hitid}};
#my@rmsdata=%rmsdata;
#print"matches : @rmsdata\n";

$hit{match}=\%rmsdata;
$hitlist[$i]=\%hit;
}
return @hitlist;
}

sub makermsfile
{
@hitlist=@{$_[0]};
my @rmspairs;


my $hitlistcount=@hitlist;
my $calc=($hitlistcount*($hitlistcount+1))/2;
print "$calc, Calculations this make take some time...\n";
my $calccount;

my $filename=$resultsdir.$runident.$rmsatomtypes."\.smr";
open(MATRIXFILE,">$filename"); # open new sequence file


print MATRIXFILE  "<XML>\n<smh type='rms'>\n";

for (my $i=0; $i<$hitlistcount; $i++)
	{
	if ($i%60==0)
		{
		$complete=(100/$calc)*$calccount;
		$complete=int($complete);
		print"\n$complete \% complete\n";
		}
	print ".";
	
	for (my $j=$i; $j<$hitlistcount; $j++)
		{
		$calccount++;
		my %hit1=%{$hitlist[$i]};
		my %hit2=%{$hitlist[$j]};
		my $ID1=$hit1{file};
		my $ID2=$hit2{file};
		#print "found hits $ID1,$ID2\n";
		my $rms=runlsqman($ID1,$ID2,$rmsatomtypes);
  	      	print "RMS:$rms $calccount of $calc\n";
		print MATRIXFILE "<RMS ID=".$ID1." ID=".$ID2 ."> ".$rms." <\\RMS>\n";
		}
}
print MATRIXFILE "<\\smh><\\XML>\n";
close MATRIXFILE;

}

1;