use clanssubs;
use File::Copy;
sub phiblastall
{
my @hitlist=@{$_[0]};
my $range=$_[1];

my $hitcount=@hitlist;

for ($i=0;$i<$hitcount;$i++)
  {
  my %hit=%{$hitlist[$i]};
  $sequence=getseq(\%hit,$range);
   
  
  }
}

sub phiblast
{
my %hit=%{$_[0]};
my $range=$_[1];#
my $Newquery=$_[2];
my $ftpdir=$_[3];
my $blastdir=$_[4];
my $seqdir=$_[5];
my $blastdb=$_[6];
my $filename=$_[7];

my $hitid=$hit{file};
my $blastout="$blastdir$hitid.phiblast";
my $basicalignment="$blastdir$hitid.align";
my $motiffile=phiblastmotif($hitid,$ftpdir,$seqdir);

if (!open(PHIFILE,$blastout))
{
print"starting blast...\n`blastpgp -i $filename  -k $motiffile  -p seedp -j 1 -o $blastout -d $blastdb`\n";
`blastpgp -i $filename  -k $motiffile  -p seedp -j 1 -o $blastout -d $blastdb`;
}
if (!open(ALIGN,$basicalignment))
{
`/home/dwoodhead/mview-1.49/bin/mview -in blast $blastout -out plain > $basicalignment`;   
}
}

sub psiblast
{
my %hit=%{$_[0]};
my $range=$_[1];
my $hitid=$hit{file};
my $blastout="/home/dwoodhead/wrkdir/blast/$hitid.psiblast";

`blastpgp -i $filename -j 2 -d $blastdb -o $blastout`
#blastpgp -i queryfilename -B alignmentfilename -j 2 -d databasefilename
#`blastpgp -i $filename  -k $motiffile  -p patseedp -o $blastout -d $blastdb`;

}

sub psipred
{
my %hit=%{$_[0]};
my $range=$_[1];
my $blastdir=$_[2];
my $seqdir=$_[3];
my $filename=$_[4];
my $hitid=$hit{file};
my $pdbid=$hit{pdbid};
$psipredout="/home/dwoodhead/$pdbid.horiz";
my $psipredfinal="$blastdir$pdbid.psipred";
print "psipred filename $filename\n";
print "psipredout filename $psipredout\n";
open (FILE,$psipredout);
@psipredtest=<FILE>;
$psipredtest=@psipredtest;
print "PSI:@psipredout\n";
print "predtest=$psipredtest\n";
if (open (BLAST,"$psipredout")&& $psipredtest!=0)
{
copy($psipredout, $psipredfinal)or die
}
else
{
print "`/home/dwoodhead/psipred/runpsipred $filename $psipredout`\n";
my @result=`/home/dwoodhead/psipred/runpsipred $filename $psipredout`;
print "Result:@result\n";
copy($psipredout, $psipredfinal)or die
#blastpgp -i queryfilename -B alignmentfilename -j 2 -d databasefilename
#`blastpgp -i $filename  -k $motiffile  -p patseedp -o $blastout -d $blastdb`;
}
}

sub runsable
{
my %hit=%{$_[0]};
my $range=$_[1];
my $seqdir=$_[3];
my $sabledir=$_[4];
my $filename=$_[5];

my $hitid=$hit{file};
$filename=~m/(.*\/)([^\/]*)\.fas/ ;
$dir=$1;
$ident=$2;
print "calculating solvent accessability...\n`/home/dwoodhead/sable_distr/run.sable $filename $ident $sabledir`\n
";
`/home/dwoodhead/sable_distr/run.sable $filename $ident $sabledir`
}

sub parsesable
{
#print "processing solvent accessability...\n";
my %hit=%{$_[0]};
my $hitid=$hit{file};
my $sabledir=[2];
my $file=$sabledir.$hitid."_RES";
$hitid=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $unknowncnt;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
my $surround=$_[1];

#print "get:$pdbid,$chain,$motifstart,$motifend\n"; 
my $openpdb="/home/dwoodhead/wrkdir/ftpdir/$pdbid.pdb";
my $seqstart=getstartresidue($openpdb,$chain);
my $offset=$motifstart - $seqstart-$surround;
#print "$offset = $motifstart - $seqstart -1 - $surround";

if ($offset<0)
{
$offset=0;
$unknowncnt=0-$offset;
}
my $length=$motifend-$motifstart+1+($surround*2);
my @outarray;

#print "opening solvent accessability file:\n $file\n";
open(AACFILE, $file);

@accfile=<AACFILE>;
close AACFILE;
$accfilecnt=@accfile;
my $seqline='X' x $unknowncnt;
my $accline='X' x $unknowncnt;
my $confline='X' x $unknowncnt;
my $getseq=0;
my $getacc=0;
my $getconf=0;

 for (my $i=0; $i<$accfilecnt; $i++)
{
	$line=$accfile[$i];
	if ($line=~m/^>.*/)	
	{
	#print "> found $line\n";
	$getseq=1;
	}
	elsif ($getseq)
	{
	#print "getseq: $line\n";
	$line=~s/\s//g;
	$seqline=$seqline.$line;
	$getseq=0;
	$getacc=1;
	}
	elsif ($getacc)
	{
	#print "getacc: $line\n";
	$line=~s/\s//g;
	$accline=$accline.$line;
	$getacc=0;
	$getconf=1;
	}
	elsif ($getconf)
	{
	#print "getconf: $line\n";
	$line=~s/\s//g;
	$confline=$confline.$line;
	$getsconf=0;
	}
}
#print "SEQ:$seqline\n";
$seqline=substr($seqline,$offset,$length);
#print "SEQ:$seqline\n";
	$seqline=~s/\n//;
	my @checkseqline=split(//,$seqline);
	my $chkseqline=@checkseqline;
	my $checkseq= $length-$chkseqline;
	#print "CHECK: $check = $length - $chkline\n";
	if ($checkseq>0)
	{
	my $endpad='X' x $checkseq;
	$seqline=$seqline.$endpad;
	}

$accline=substr($accline,$offset,$length);
#print"ACC:$accline\n";
	$accline=~s/\n//;
	my @checkaccline=split(//,$accline);
	my $chkaccline=@checkaccline;
	my $checkacc= $length-$chkaccline;
	#print "CHECK: $check = $length - $chkline\n";
	if ($checkacc>0)
	{
	my $endpad='0' x $checkacc;
	$accline=$accline.$endpad;
	}
$confline=substr($confline,$offset,$length);
#print "CONF:$confline\n";
$confline=~s/\n//;
	my @checkconfline=split(//,$confline);
	my $chkconfline=@checkconfline;
	my $checkconf= $length-$chkconfline;
	#print "CHECK: $check = $length - $chkline\n";
	if ($checkconf>0)
	{
	my $endpad='0' x $checkconf;
	$confline=$confline.$endpad;
	}




my @seqline=split(//,$seqline);
my @accline=split(//,$accline);
my @confline=split(//,$confline);

my $seqcnt=@seqline;

for (my $i=0; $i<$seqcnt; $i++)
	{
#	print "SABLE:POS:$i-$surround,RES:$seqline[$i],ACC:$accline[$i],CONF:$confline[$i]\n";
	my %res=('pos',$i-$surround,'res',$seqline[$i],'acc',$accline[$i],'conf',$confline[$i]);
	push (@outarray, \%res);
	}
return @outarray;
}

sub solacc
{
my %hit=%{$_[0]};
my $range=$_[1];
my $hitid=$hit{file};
my $out="/home/dwoodhead/wrkdir/sequences/$hitid.acc";
my $filename=saveseqfile(\%hit,$range,$Newquery);
#print "calculating solvent accessability...\n`/home/dwoodhead/sspro4/bin/predict_acc.sh $filename $out`\n";
`/home/dwoodhead/sspro4/bin/predict_acc.sh $filename $out`
}

sub calcacc
{
print "processing solvent accessability...\n";
my %hit=%{$_[0]};
my $hitid=$hit{file};
my $file="/home/dwoodhead/wrkdir/sequences/$hitid.acc";
$hitid=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
#print "get:$pdbid,$chain,$motifstart,$motifend\n"; 
my $openpdb="/home/dwoodhead/wrkdir/ftpdir/$pdbid.pdb";
my $seqstart=getstartresidue($openpdb,$chain);
my $offset=$motifstart - $seqstart-1-$surround;
my $length=$motifend-$motifstart+1+($surround*2);

open(AACFILE, $file);

@aacfile=<AACFILE>;
close AACFILE;

$line=$aacfile[0];
$line2=$aacfile[1];
$line3=$aacfile[2];

$motif2=substr($line2,$offset,$length);
#print "ACC2:$motif2\n";
$motif3=substr($line3,$offset,$length);
print"$motif2\n$motif3\n";
#print "ACC3:$motif3\n";

return $motif3;
}

sub acccal
{
my $motif3=$_[0];
my @motif3=split(//,$motif3);

my $accscore=0;
for $each(@motif3)
{
$accscore=$accscore+accscore($each);
}
#print "$accscore\n";
return $accscore;
}

sub accscore
{
my $letter=$_[0];

if ($letter eq 'b')
{return 0}
if ($letter eq 'e')
{return 1}
}

sub accscore
{
my $letter=$_[0];

if ($letter eq 'b')
{return 0}
if ($letter eq 'e')
{return 1}
}

sub altaccscore
{
my $letter=$_[0];

if ($letter eq 'b')
{return -1}
if ($letter eq 'e')
{return 1}
}

sub getnextconserved
{
my $id=$_[0];
my $blastdir=$_[2];
my $file="$blastdir$id.align";
my $surround=$_[1];
#print "$id\n";
$id=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
#print "getcon:$pdbid,$chain,$motifstart,$motifend\n"; 
my $openpdb="/home/dwoodhead/wrkdir/ftpdir/$pdbid.pdb";
my $seqstart=getstartresidue($openpdb,$chain);

my $offset=$motifend - $seqstart +1;
#print "$motifend - $surround - $seqstart -1= $offset\n";
#my $length=$motifend-$motifstart+($surround)+1;
#print "$motifend - $motifstart + ($surround)\n";
my @conresarray;


open (ALIGNMENT ,"<$file");
my @alignment=<ALIGNMENT>;
#print "ALIGN: @alignment";
close ALIGNMENT;
my $linecnt=@alignment;
for my $line1(@colinfo)
{
#print "$line1\n";
}

my @alignmentarray;
for (my $i=0; $i<$linecnt; $i++)
	{
	my $line=@alignment[$i];
	my @line=split (//,$line);
	my $length=@line;
	$line=~s/\S*\s*//;
#	print "New:$line";
	#print "$offset, $length\n";
	my $newline=substr($line,$offset,$length);
#	print "substr:$newline\n";
	push(@alignmentarray,$newline);
	}


my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;

while($checkres<$firstlngth)
{
#print "conserved res not found, findconservedres:$findconservedres outarraycount:$outarraycnt\n";
if($firstline[$checkres]eq"D"||$firstline[$checkres]eq"E")
{ 
my $matcharraycnt=$outarraycnt;
my $conresnum=$checkres;
#print "found($checkres) res $conresnum:$firstline[$checkres]\n ";
my $matchcntD=0;
my $matchcntE=0;
	for my $each(@alignmentarray)
		{
		
		#print "012345678901234567890123456789\n";
		#print "$each\n";
		my @each=split(//,$each);
		#print "$checkres, $firstline[$checkres]:$each[$checkres],\n\n";
		if ($each[$checkres]eq "D")
			{
			++$matchcntD;
			}
		if($each[$checkres]eq"E")
			{
			++$matchcntE;
			}
		if($each[$checkres] eq "-")
			{
			--$matcharraycnt;
			}
		else
			{}
		}
	my $percentconD=($matchcntD/$matcharraycnt)*100;
	my $percentconE=($matchcntE/$matcharraycnt)*100;
	my $totpercentcon=$percentconD+$percentconE;
	#print"$totpercentcon=$percentconD+$percentconE\n";
	
	if ($totpercentcon>50)
		{
#		print "$conresnum:$totpercentcon % conserved: E:$percentconE D:$percentconD\n";
		my %conservedres=('dist',$conresnum,'percentD',$percentconD,'percentE',$percentconE);
		push (@conresarray, \%conservedres);
		#$conservedresfound=1
		}
	else
		{}
}
else
{
}
++$checkres;
}
return @conresarray;
}

sub testconserved
{
my $top=$_[1];
my $bottom=$_[2];
my @conserved=@{$_[0]};
my $concnt=@conserved;
my $conserved=0;
	for (my $i=0; $concnt>$i; $i++)
	{
	my %con=%{$conserved[$i]};
	#my @con=%con;
	#print "ConData:@con\n";
	my $dist=$con{dist};
	if
	($dist<=$top && $dist>= $bottom)
		{	
		$conserved++;
		}
	else	{}
	}
if ($conserved!=0)
{return "True";}
else
{return "False";}
}

sub testconserved2
{
my $top=$_[1];
my $bottom=$_[2];
my @conserved=@{$_[0]};
my $concnt=@conserved;
my $conserved=0;
	for (my $i=0; $concnt>$i; $i++)
	{
	my %con=%{$conserved[$i]};
	#my @con=%con;
	#print "ConData:@con\n";
	my $dist=$con{dist};
	if
	($dist<=$top && $dist>= $bottom)
		{	
		$conserved++;
		}
	else	{}
	}
if ($conserved!=0)
{return "1";}
else
{return "-1";}
}

sub getcolinfo
{
my $id=$_[0];
my $blastdir=$_[2];
my $file="$blastdir$id.align";
my $surround=$_[1];
#print "$id\n";
$id=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
my $unknowncnt;
#print "get:$pdbid,$chain,$motifstart,$motifend\n"; 
my $openpdb="/home/dwoodhead/wrkdir/ftpdir/$pdbid.pdb";
my $seqstart=getstartresidue($openpdb,$chain);

my $offset=$motifstart -$surround - $seqstart;

#print "$motifstart - $surround - $seqstart = $offset\n";
my $unknowncnt=0-$offset;
if ($offset<0)
{$offset=0;}
my $length=$motifend-$motifstart+(2*$surround)+1;
#print "$length = $motifend - $motifstart + (2x$surround)\n";

my @outarray;

#print "file:$file\n";
open (ALIGNMENT ,"<$file");
my @alignment=<ALIGNMENT>;
#print "alignment:@alignment\n";
close ALIGNMENT;
my $linecnt=@alignment;

for (my $i=0; $i<$linecnt; $i++)
	{
	$line=@alignment[$i];
	$line=~s/[\w,|,.]*\s*//;
	#print "New:$line";
	#print "$offset, $length\n";
	$newline='X' x  $unknowncnt;
	$newline=$newline.substr($line,$offset,$length);
	#print "newline:\"$newline\"\n";
	$newline=~s/\n//;
	my @checkline=split(//,$newline);
	my $chkline=@checkline;
	my $check= $length-$chkline;
	#print "CHECK: $check = $length - $chkline\n";
	if ($check>0)
	{
	my $endpad='X' x $check;
	$newline=$newline.$endpad;
	}
	#print "substr:$newline\n";
	push(@outarray,$newline);
	#$i++;
	}
#print "COL outarray: @outarray\n";
return @outarray;
}

sub processcoldata
{
#print "processcldata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];
my $ss=$_[2];

while($checkres<$firstlngth)
{
my $sizearraycnt=$outarraycnt;
my $sizecntS=0;
my $sizecntL=0;
my $sizecntU=0;
my $totalsize=0;


my $typearraycnt=$outarraycnt;
my $tcntPl=0;
my $tcntNe=0;
my $tcntP=0;
my $tcntN=0;
my $tcntG=0;
my $tcntU=0;
my %aatype;
my $meanarraycnt=$outarraycnt;

	for (my $i=0;$i<$outarraycnt;$i++)
		{
		my $each=@alignmentarray[$i];
		#print "012345678901234567890123456789\n";
		#print "$each\n";
		my @each=split(//,$each);
		my $aaid=$each[checkres];
		++$aatype{'$aaid'};
		
		my $eachsize=sortressize($each[$checkres],$ss);
		my $eachtype=sortrestype($each[$checkres]);
		
		my $weight=aasize($each[$checkres]);
		if ($weight)
		{
		$totalsize=$totalsize+$weight;
		}
		else
			{
			--$meanarraycnt;
			#print "*$eachsize*\n";
			}

		#print "$checkres,$each[$checkres]:$eachsize:$eachtype,$totalsize\n ";
		if ($eachsize eq "S")
			{
			++$sizecntS;
			}
		elsif($eachsize eq "L")
			{
			++$sizecntL;
			}
		else
			{
			++$sizecntU;
			--$sizearraycnt;
			#print "*$eachsize*\n";
			}

		if ($eachtype eq "+")
			{
			++$tcntPl;
			}
		elsif($eachtype eq "-")
			{
			++$tcntNe;
			}
		elsif($eachtype eq "P")
			{
			++$tcntP;
			}
		elsif ($eachtype eq "N")
			{
			++$tcntN;
			}
		elsif($eachtype eq "G")
			{
			++$tcntG;
			}
		else
			{
			++$tcntU;
			--$typearraycnt;
			#print "*$eachtype*\n";
			}
		}
#	print "$sizecntS/$sizearraycnt*100;\n";
	my $percentconS;
	my $percentconL;
	my $percentconU;
	my $perPl;
	my $perNe;
	my $perP;
	my $perN;
	my $perG;
	my $perU;
	if ($sizearraycnt)
	{
	$percentconS=$sizecntS/$sizearraycnt*100;
	$percentconL=$sizecntL/$sizearraycnt*100;
	$percentconU=$sizecntU/$sizearraycnt*100;
	}
	else
	{
	$percentconS=0;
	$percentconL=0;
	$percentconU=0;
	
	}
	
	if ($meanarraycnt)
	{
	$meansize=$totalsize/$meanarraycnt;
	#print "Meansize:$totasize\/$sizearraycnt\=$meansize\n";
	}
	else
	{
	$meansize=0;
	#print "Meansize:$totasize\/$sizearraycnt\=$meansize\n";
	}
	
	if ($typearraycnt)
	{
	$perPl=$tcntPl/$typearraycnt*100;
	$perNe=$tcntNe/$typearraycnt*100;
	$perP=$tcntP/$typearraycnt*100;
	$perN=$tcntN/$typearraycnt*100;
	$perG=$tcntG/$typearraycnt*100;
	$perU=$tcntU/$typearraycnt*100;
	
	}
	else
	{
	$perPl=0;
	$perNe=0;
	$perP=0;
	$perN=0;
	$perG=0;
	$perU=0;
	}
	#print"$percentcon=$matchcnt/$outarraycnt*100\n";
	my $conresnum=$checkres-$surround;
	#print "POSITION: $conresnum\n ";
	#print "AASIZE S:$percentconS L:$percentconL\n"; 
	#print " U:$percentconU ";
	#print "AATYPE +:$perPl -:$perNe P:$perP N:$perN G:$perG\n";
	#print " U:$perU\n";
	#print "meansize:$meansize\n";
	my %conservedres=('pos',$conresnum,'res',$each[$checkres],'perL', $percentconL,'perS',$percentconS, 'perPl',$perPl,'perNe',$perNe,'perP',$perP,'perN',$perN,'perG',$perG,'meansize',$meansize);
	if ($typearraycnt)
	{
	my @aatype=%aatype;
	my $aatypecnt;
	for (my $i=0;$i<$aatypecnt;$i=$i+2)
		{
		my $id=$aatype[$i];
		my $aacnt=$aatype[$i+1];
		$aacnt=($aacnt/$outarraycnt)*100;
		print "$id:$aacnt\n";
		$conservedres{$id}=$aacnt;
		}
	}
	push (@conresarray, \%conservedres);
	#$conservedresfound=1

	

++$checkres;
}
return @conresarray;

}

sub processcoldatah
{
#print "processcldata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];
my $ss=$_[2];

while($checkres<$firstlngth)
{
my $sizearraycnt=$outarraycnt;
my $sizecntPI=0;
my $sizecntPO=0;
my $sizecntU=0;
my $totalhyd=0;
my $meanarraycnt=$outarraycnt;

	for (my $i=0;$i<$outarraycnt;$i++)
		{
		my $each=@alignmentarray[$i];
		#print "012345678901234567890123456789\n";
		#print "$each\n";
		my @each=split(//,$each);
		
	
		my $aahyd=aahyd($each[$checkres]);
		if ($aahyd)
		{
		$totalhyd=$totalhyd+$aahyd;
		}
		else
			{
			--$meanarraycnt;
			#print "*$eachsize*\n";
			}

		my $eachsize=sorthyd($each[$checkres],$ss);
		
		#print "$checkres,$each[$checkres]:$eachsize\n ";
		if ($eachsize eq "PHIL")
			{
			++$sizecntPI;
			}
		elsif($eachsize eq "PHOB")
			{
			++$sizecntPO;
			}
		else
			{
			++$sizecntU;
			--$sizearraycnt;
			#print "*$eachsize*\n";
			}

		
		}
#	print "$sizecntS/$sizearraycnt*100;\n";
	my $percentconPI;
	my $percentconPO;
	my $percentconU;
	
	if ($sizearraycnt)
	{
	#print "**sizearraycnt$sizearraycnt\n";
	$percentconPI=$sizecntPI/$sizearraycnt*100;
	$percentconPO=$sizecntPO/$sizearraycnt*100;
	$percentconU=$sizecntU/$sizearraycnt*100;
	}
	else
	{
	#print "**sizearraycnt FALSE\n";
	$percentconPI=0;
	$percentconPO=0;
	$percentconU=0;
	}
	if ($meanarraycnt)
	{
	$meanhyd=$totalhyd/$meanarraycnt;
	#print "Meansize:$totasize\/$sizearraycnt\=$meansize\n";
	}
	else
	{
	$meanhyd=0;
	#print "Meansize:$totasize\/$sizearraycnt\=$meansize\n";
	}
	
	#print"$percentcon=$matchcnt/$outarraycnt*100\n";
	my $conresnum=$checkres-$surround;
	#print "meanhyd:$meanhyd\n";
	#print "POSITION: $conresnum\n ";
	#print "AAHYD PI:$percentconPI PO:$percentconPO\n"; 
	#print " U:$percentconU ";
	#print "AATYPE +:$perPl -:$perNe P:$perP N:$perN G:$perG\n";
	#print " U:$perU\n";
	my %conservedres=('pos',$conresnum,'res',$each[$checkres],'perPI', $percentconPI,'perPO',$percentconPO,'meanhyd',$meanhyd,);
	push (@conresarray, \%conservedres);
	#$conservedresfound=1

	

++$checkres;
}
return @conresarray;

}

sub getaadata
{
#print "getabsolutedata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];

while($checkres<$firstlngth)
{
my $each=@alignmentarray[1];
my @each=split(//,$each);
		

my $eachtype=$each[$checkres];

my $conresnum=$checkres-$surround;
#print "ABSRES : $conresnum, $eachsize, $eachtype\n";
my %conservedres=('pos',$conresnum,'AA',$eachtype);
push (@conresarray, \%conservedres);
	#$conservedresfound=1
++$checkres;
}
return @conresarray;
}

sub getaasize
{
#print "getabsolutedata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];

while($checkres<$firstlngth)
{
my $each=@alignmentarray[1];
my @each=split(//,$each);
		

my $eachtype=aasize($each[$checkres]);

my $conresnum=$checkres-$surround;
#print "SIZE : $conresnum, $eachtype\n";
my %conservedres=('pos',$conresnum,'size',$eachtype);
push (@conresarray, \%conservedres);
	#$conservedresfound=1
++$checkres;
}
return @conresarray;
}

sub getaahyd
{
#print "getabsolutedata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];

while($checkres<$firstlngth)
{
my $each=@alignmentarray[1];
my @each=split(//,$each);
		

my $eachtype=aahyd($each[$checkres]);

my $conresnum=$checkres-$surround;
#print "SIZE : $conresnum, $eachtype\n";
my %conservedres=('pos',$conresnum,'aahyd',$eachtype);
push (@conresarray, \%conservedres);
	#$conservedresfound=1
++$checkres;
}
return @conresarray;
}
sub getabsolutedata
{
#print "getabsolutedata...\n";
my @alignmentarray=@{$_[0]};
#print ":@alignmentarray\n";
my $checkres=0;
my $firstline=$alignmentarray[0];
my $outarraycnt=@alignmentarray;
#print "Firstline: $firstline\n";
my @firstline=split(//,$firstline);
my $firstlngth=@firstline;
my @conresarray;
my $surround=$_[1];

while($checkres<$firstlngth)
{
my $each=@alignmentarray[1];
my @each=split(//,$each);
		
my $eachsize=sortressize($each[$checkres]);
my $eachtype=sortrestype($each[$checkres]);

my $conresnum=$checkres-$surround;
#print "ABSRES : $conresnum, $eachsize, $eachtype\n";
my %conservedres=('pos',$conresnum,'size',$eachsize,'type',$eachtype);
push (@conresarray, \%conservedres);
	#$conservedresfound=1
++$checkres;
}
return @conresarray;
}

sub alignsecondary
{
my $id=$_[0];
my $blastdir=$_[1];
my $psipreddir=$_[2];
my $file="$psipreddir$id.horiz";
open (FILE , $file);
my @psipredfile=<FILE>;
close FILE;

my $secondaypred;
my $seq;

for $line(@psipredfile)
{
if ($line=~m/Pred:\s(\w*)/)
{
$secondarypred=$secondarypred.$1
}
if ($line=~m/AA:\s(\w*)/)
{
$seq=$seq.$1
}
else
{}

}
my $fileout="$blastdir$id.pred";

my $out=$secondarypred."\n".$seq;
#print $out."\n";
open (FILEOUT, ">$fileout");
print FILEOUT $out;
close FILEOUT;
}

sub secondarydata
{
my $id=$_[0];
my $psipreddir=$_[1];
my $file="$psipreddir$id.pred";

#print "$id\n";
$id=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
#print "get:$pdbid,$chain,$motifstart,$motifend\n"; 
my $openpdb="/home/dwoodhead/wrkdir/ftpdir/$pdbid.pdb";
my $seqstart=getstartresidue($openpdb,$chain);

my $offsets=$motifstart - $seqstart-1;
#print "Motifstart:$motifstart - seqstart:$seqstart +1 = offset:$offset\n";
my $unknowncnt=0-$offsets;
if ($offsets<0)
{$offsets=0;}

my $length=$motifend-$motifstart+1;
#print "motifend:$motifend - motifstart:$motifstart +1 = length:$length\n";

my $offsete=$offsets+$length-$unknowncnt;


open(PREDFILE, $file);

@predfile=<PREDFILE>;
close PREDFILE;

my $line=$predfile[0];
my $line2=$predfile[1];

my $testlineup=substr($line,0,$offsets);
#print "utestline:$testlineup\n";
$testlineup=~m/[C]([H,E]*)([C]*)$/;
#print "udistline:$2\n";
my @uparray=split (//, $2);
my $updist=@uparray;

$updist=$updist;
#print "";
@usslength=split (//, $1);
$usslength=@usslength;
$uptype=$usslength[$uplenth-1];
if ($uptype eq "")
{$uptype="-";}
#print "usslengthline:$1\nusslength:$usslength\n$uptype\n\n";

my $testlinemotif=substr($line,$offsets,$length);
#print "mtestline:$testlinemotif\n\n";
my $ssscore=0;
my @mtest=split(//,$testlinemotif);
for my $each(@mtest)
{
if ($each eq H||$each eq E)
{
--$ssscore;
}
else{}
}
#print "ss Score:$ssscore\n";

my $testlinedown=substr($line,$offsete);
#print "dtestline:$testlinedown";
$testlinedown=~m/^([C]*)([H,E]*)[C*]/;

#print "ddistline:$1\n";
my @downarray=split (//, $1);
my $downdist=@downarray;

#print "Down dist= $downdist\n";

#print "dsslenghtline:$2\n";
@dsslength=split (//, $2);
$dsslength=@dsslength;
$downtype=$dsslength[0];

if ($downtype eq "")
{$downtype="-";}
#print "dssline:$2\ndsslength:$dsslength\n$downtype\n\n";


my $getstart=$offsets-$updist-$usslength-1;
my $getlength=$updist+$usslength+5+$dsslength+$downdist+2;
my $getline1=substr($line,$getstart,$getlength);
my $getline2=substr($line2,$getstart,$getlength);

my %outdata=('updist',$updist,'uplength',$usslength,'uptype',$uptype,'ssscore',$ssscore,'downdist',$downdist,'downlength',$dsslength,'downtype',$downtype,'ssseq1',$getline1,'ssseq2',$getline2)
}

sub sortressize
{

my $residue=$_[0];
$residue=uc($residue);
my $ss=$_[1];
if ($ss==1)
{
if($residue eq G){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==2)
{
if($residue eq G||$residue eq A){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==3)
{
if($residue eq G||$residue eq A||$residue eq S){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==4)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==5)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==6)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==7)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==8)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq N||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==9)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W||$residue eq D) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==10)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D){return "S";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==11)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K){return "S";}
elsif($residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==12)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E){return "S";}
elsif($residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==13)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E||$residue eq M){return "S";}
elsif($residue eq H||$residue eq F||$residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==14)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H){return "S";}
elsif($residue eq F||$residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==15)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F){return "S";}
elsif($residue eq R||$residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==16)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R){return "S";}
elsif($residue eq Y||$residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==17)
{
if($residue eq G||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D||$residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq R||$residue eq Y){return "S";}
elsif($residue eq W) {return "L";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
}
sub sorthyd
{

my $residue=$_[0];
$residue=uc($residue);
my $ss=$_[1];
if ($ss==1)
{
if($residue eq R){return "PHIL";}
elsif($residue eq Q||$residue eq K||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq G||$residue eq Y||$residue eq W||$residue eq A||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==2)
{
if($residue eq R||$residue eq K){return "PHIL";}
elsif($residue eq Q||$residue eq G||$residue eq E||$residue eq M||$residue eq H||$residue eq F||$residue eq A||$residue eq Y||$residue eq W||$residue eq S||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq N||$residue eq D) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==3)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq H||$residue eq F||$residue eq G||$residue eq Y||$residue eq W||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq S) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==4)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq Y||$residue eq W||$residue eq P||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq S) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==5)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq Y||$residue eq W||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq S) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==6)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq W||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq S) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==7)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L||$residue eq S) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==8)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq V||$residue eq T||$residue eq C||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==9)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq G||$residue eq V||$residue eq C||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==10)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G){return "PHIL";}
elsif($residue eq A||$residue eq M||$residue eq F||$residue eq V||$residue eq C||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==11)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A){return "PHIL";}
elsif($residue eq M||$residue eq F||$residue eq V||$residue eq C||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==12)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A||$residue eq M){return "PHIL";}
elsif($residue eq F||$residue eq V||$residue eq C||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==13)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A||$residue eq M||$residue eq C){return "PHIL";}
elsif($residue eq F||$residue eq V||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==14)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A||$residue eq M||$residue eq C||$residue eq F){return "PHIL";}
elsif($residue eq V||$residue eq I||$residue eq L) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==15)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A||$residue eq M||$residue eq C||$residue eq F||$residue eq L){return "PHIL";}
elsif($residue eq V||$residue eq I) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
if ($ss==16)
{
if($residue eq R||$residue eq K||$residue eq N||$residue eq D||$residue eq Q||$residue eq E||$residue eq H||$residue eq P||$residue eq Y||$residue eq W||$residue eq S||$residue eq T||$residue eq G||$residue eq A||$residue eq M||$residue eq C||$residue eq F||$residue eq V||$residue eq L){return "PHIL";}
elsif($residue eq I) {return "PHOB";}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") {return "X";}
else{print "**$residue**\n";}
}
}
sub sortrestype
{
my $residue=$_[0];
$residue=uc($residue);
if($residue eq 'R'||$residue eq 'H'||$residue eq 'K') 
{
return "+";
}
elsif($residue eq 'D' ||$residue eq 'E')
{
return "-";
}
elsif($residue eq 'S'||$residue eq 'T'||$residue eq 'N'||$residue eq 'Q'||$residue eq 'C')
{
return "P";
}
elsif($residue eq A||$residue eq I||$residue eq L||$residue eq M||$residue eq F||$residue eq W||$residue eq Y||$residue eq V||$residue eq P)
{
return "N";
}
elsif($residue eq G)
{
return "G";
}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") 
{
return "X";
}
else
{
#print "**$residue**\n";
}
}
sub aasize
{
my $residue=$_[0];
$residue=uc($residue);

if($residue eq 'G') 
{
return "75";
}
elsif($residue eq 'A') 
{
return "89";
}
elsif($residue eq 'S') 
{
return "105";
}
elsif($residue eq 'P') 
{
return "115";
}
elsif($residue eq 'V')
{
return "117";
}
elsif($residue eq 'T')
{
return "119";



}
elsif($residue eq 'C')
{
return "121";
}
elsif($residue eq 'L'||$residue eq 'I')
{
return "131";
}
elsif($residue eq 'N')
{
return "132";
}
elsif($residue eq 'D')
{
return "133";
}
elsif($residue eq 'Q'||$residue eq 'K')
{
return "146";
}
elsif($residue eq 'E')
{
return "147";
}
elsif($residue eq 'M')
{
return "149";
}
elsif($residue eq 'H')
{
return "155";
}
elsif($residue eq 'F')
{
return "165";
}
elsif($residue eq 'R')
{
return "174";
}
elsif($residue eq 'Y')
{
return "181";
}
elsif($residue eq 'W')
{
return "204";
}

elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") 
{
return "0";
}
else
{
#print "**$residue**\n";
}
}
sub aahyd
{
my $residue=$_[0];
$residue=uc($residue);

if($residue eq 'G') 
{
return "-0.4";
}
elsif($residue eq 'A') 
{
return "1.8";
}
elsif($residue eq 'S') 
{
return "-0.8";
}
elsif($residue eq 'P') 
{
return "-1.6";
}
elsif($residue eq 'V')
{
return "4.2";
}
elsif($residue eq 'T')
{
return "-0.7";
}
elsif($residue eq 'C')
{
return "2.5";
}
elsif($residue eq 'L')
{
return "3.8";
}
elsif($residue eq 'I')
{
return "4.5";
}
elsif($residue eq 'N'||$residue eq 'D'||$residue eq 'E'||$residue eq 'Q')
{
return "-3.5";
}
elsif($residue eq 'K')
{
return "-3.9";
}
elsif($residue eq 'M')
{
return "1.9";
}
elsif($residue eq 'H')
{
return "-3.2";
}
elsif($residue eq 'F')
{
return "2.8";
}
elsif($residue eq 'R')
{
return "-4.5";
}
elsif($residue eq 'Y')
{
return "-1.3";
}
elsif($residue eq 'W')
{
return "-0.9";
}

elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") 
{
return "0";
}
else
{
#print "**$residue**\n";
}
}
sub saveseqfile
{
my %hit=%{$_[0]};
my $range=$_[1];
my $ftpdir=$_[2];
my $seqdir=$_[3];
my $file;

my $sequence=getseq(\%hit,$range,,$ftpdir);
#print $sequence;
my $hitid=$hit{file};
#print "making sequence file for $hitid\n";
if ($range)
{$file="$seqdir$hitid[$range].fas"}
else
{$file="$seqdir$hitid.fas"}

#print "creating file :$file\n";
open (MYSEQFILE,">$file");
print MYSEQFILE "\>$hitid\n";
print MYSEQFILE "$sequence\n";
close MYSEQFILE;
return $file;
}

sub phiblastmotif
{
my $id=$_[0];
my $ftpdir=$_[1];
my $seqdir=$_[2];
$id=~m/(.{4})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifend=$4;
#print "get:$pdbid,$chain,$motifstart,$motifend\n";
my $file="$seqdir$id.motif";
my $openpdb="$ftpdir$pdbid.pdb";
#print "get Start residue...\n";
my $seqstart=getstartresidue($openpdb,$chain);
my $offset=$motifstart - $seqstart +1;
my $offset2=$offset+4;

#print "creating file :$file\n";
#print "$id: seqstart:$seqstart, motif:$motifstart\($offset\)-$motifend\($offset2\)\n"; 
open (MYMOTIFFILE,">$file");
print MYMOTIFFILE "ID  $id DXDXDG motif\n";
print MYMOTIFFILE "PA  D-x-[DNS]-x-[DNS]\n";
print MYMOTIFFILE "HI ($offset $offset2)\n";
close MYMOTIFFILE;
return $file;
}

sub checkphiblast
{
my %hit=%{$_[0]};
my $hitid=$hit{file};
my $blastout="/home/dwoodhead/wrkdir/blast/$hitid.phiblast";
my $basicalignment="/home/dwoodhead/wrkdir/blast/$hitid.align";


}

sub checkpsipred
{
my %hit=%{$_[0]};
my $range=$_[1];
my $hitid=$hit{file};
my $psipredout="/home/dwoodhead/wrkdir/blast/$hitid.psiblast";
$filename=saveseqfile(\%hit,$range);
#print "`/home/dwoodhead/psipred/runpsipred $filename`\n";

}


sub skip
{
my @families=@{$_[0]};
my $skip=$_[1];
my @outarray;
my $arraycnt=@families;

for (my $i=0;$i<$arraycnt;$i++)
{
my %hit=%{$families[$i]};
$scop=$hit{scop};
if($scop eq $skip)
{
#print "$i skipping family: $hit{scop}\n";
}
else
{
#print "$i $scop ne $skip\n";
push (@outarray, \%hit);
}
}
return @outarray;

}

sub mergefile
{
$file1=$_[0];
$file2=$_[1];
$file3=$_[2];
$outfile=$_[3];
open (OUTFILE, ">$outfile");
open (FILE1,"<$file1");
open (FILE2,"<$file2");
open (FILE3,"<$file3");
@file1=<FILE1>;
@file2=<FILE2>;
@file3=<FILE3>;
for my$line(@file1)
{print OUTFILE $line;}
print OUTFILE "\n";
for my$line(@file2)
{print OUTFILE $line;}
print OUTFILE "\n";
for my$line(@file3)
{print OUTFILE $line;}
print OUTFILE "\n";
}

sub getunique
{
my @thresholds=@{$_[0]};
my @outlist;

#print "thresholds::@thresholds\n";

for my $each(@thresholds)
	{
	my $flag=0;
	for my $every(@outlist)
	{
	if ($each == $every)
		{
		#print "$each == $every\n";
		$flag=1;
		}
	
	}
	if ($flag==0 && $each!=0)
		{
		#print "new threshold $each";
		push (@outlist,$each);
		}
	
	}
#print "threshold list:: @outlist\n";
return @outlist;
}

sub writeerrordata
{
my $file=$_[1];
my %errordata=%{$_[0]};
my $errorfile=$_[2];
my $skip=$_[3];
open (ERRORFILE,">>$errorfile");

open (FILE,">>$file");
my @errors=@{$errordata{errors}};
my $per=$errordata{pererror};
my $ernum=$errordata{numerror};
my $ercnt=@errors;
print ERRORFILE "$skip\t$per\t$ernum\t";
print DTFILE "ERRORS\n";
for ($i=0;$i<$ercnt;$i++)
{
my %hit=%{$errors[$i]};
my $dtid=$hit{dtid};
my $id=$hit{repid};
my $scop=$hit{scop};
print FILE "$dtid\t$id\t$scop\t";
print ERRORFILE "$dtid\t$id\t$scop\t";

}
print ERRORFILE "\n";

}

sub writeSVMdatafile
{

my %errordata=%{$_[0]};
my $file=$_[1];
my $errorfile=$_[2];
open (ERRORFILE,">>$errorfile");

my $miss=$errordata{miss};
my $errorest=$errordata{errorest};
my $recallest=$errordata{recallest};
my $precisoinest=$errordata{precisionest};

print ERRORFILE "$miss\t$errorest\t$recallest\t$precisionest\t";
print ERRORFILE "\n";

}
1;
