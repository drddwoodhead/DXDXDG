sub getfamilynames
{
my $startfile=$_[0];

my @families=parsexml($startfile);
my @outarray;
#my $familycnt=8;
my $familycnt=@families;
my $surround=$_[1];
for(my $i=0;$i<$familycnt;$i++)
	{
	my %family=%{$families[$i]};
	my @idarray=@{$family{pdblist}};
	my $idcnt=@idarray;
	
	for (my $j=0; $j<$idcnt; $j++)
	{	
	my %newfamily=%family;
	my $id=$idarray[$j];
	$newfamily{'repid'}=$id;
	my %hit=temphithash($id);
	my $file="/home/dwoodhead/wrkdir/ftpdir/$id.pdb";
	push (@outarray, \%newfamily);
	} 
	}
return @outarray;
}

sub processfamilydata
{

my $startfile=$_[0];

my @families=parsexml($startfile);
my @outarray;
#my $familycnt=8;
my $familycnt=@families;
my $surround=$_[1];
my $aas=$_[2];
my @aathresholdlist=@{$_[3]};
my $aah= $_[4];
my @aahthresholdlist=@{$_[5]};
my $con=$_[6];
my $sol=$_[7];
my $ss=$_[8];
my $blast=$_[9];
my $blastdir=$_[10];
my $psipreddir=$_[11];
my $sabledir=$_[12];


for(my $i=0;$i<$familycnt;$i++)
	{
	my %family=%{$families[$i]};
	my @idarray=@{$family{pdblist}};
	my $idcnt=@idarray;
	$bvalue=$family{bvalue};
	
	for (my $j=0; $j<$idcnt; $j++)
	{	
	my %newfamily=%family;
	my $id=$idarray[$j];
	$newfamily{'repid'}=$id;
	$newfamily{'bvalue'}=$bvalue;
	my %hit=temphithash($id);
	#my $file="/home/dwoodhead/wrkdir/ftpdir/$id.pdb";

	#print "\nprocessing ID:$id SCOP:".$family{scop}." files...family :".($i+1)." of $familycnt hit:".($j+1)." of $idcnt\n";
	$newfamily{'alignfile'}=$blastdir.$id.".align";	
	
	if ($blast!=0)
	{
	#print "processing blast data...\n";
	my @colinfo=getcolinfo($id,$surround, $blastdir);
	
	if ($aas eq "v")
	{
	#print "aas = v\n";
	my @thresholds=getunique(\@aathresholdlist);
	my $thresholdnum=@thresholds;	

	for(my $i=0; $i<$thresholdnum;$i++)
	{
	my $aath=$thresholds[$i];
	#print "getting aas data for threshold $aath\n";
	my @restypes=processcoldata(\@colinfo,$surround,$aath);
	$newfamily{'resdata'.$aath}=\@restypes;
	
	}
	}	
	

	#print "aas \!\=v\n";
	my @restypes=processcoldata(\@colinfo,$surround,$aas);
	$newfamily{'resdata'}=\@restypes;
	

	if ($aah eq "v")
	{
	#print "aas = v\n";
	my @thresholds=getunique(\@aahthresholdlist);
	my $thresholdnum=@thresholds;	

	for(my $i=0; $i<$thresholdnum;$i++)
	{
	my $aath=$thresholds[$i];
	my @restypes=processcoldatah(\@colinfo,$surround,$aah);
	$newfamily{'resdatah'.$aath}=\@restypes;
	}
	}	
	else
	{
	#print "aas /=v\n";
	my @restypes=processcoldatah(\@colinfo,$surround,$aah);
	$newfamily{'resdatah'}=\@restypes;
	}

	my @absolutedata=getabsolutedata(\@colinfo,$surround);
	my @aadata=getaadata(\@colinfo,$surround);
	my @aasize=getaasize(\@colinfo,$surround);
	my @aahyd=getaahyd(\@colinfo,$surround);
	
	$newfamily{'absdata'}=\@absolutedata;
	$newfamily{'AAdata'}=\@aadata;
	$newfamily{'sizedata'}=\@aasize;
	$newfamily{'hyddata'}=\@aahyd;
	#print "RESDATA:@restypes\n";
	}

	if ($ss!=0)
	{
	#print "processing psipred files...\n";
	alignsecondary($id,,$psipreddir);
	$newfamily{'SSfile'}=$blastdir.$id.".pred";
	my %ssdata=secondarydata($id,$psipreddir);
	$newfamily{'SSdata'}=\%ssdata;
	#print "SSDATA:%ssdata\n";
	}
	
	if ($con!=0)
	{
	#print "finding conserved downstream D/E...\n";
	my @conserved=getnextconserved($id,99,$blastdir);
	$newfamily{'conserved'}=\@conserved;
	#print "CONDE:@conserved\n";
	}
	if ($sol!=0)
	{
	#print "processing solvent accesability files...\n";
	my @accarray=parsesable(\%hit,$surround,$sabledir);
	my $test=@aacarray;
	#print "**AACARRAY::$test\n";
	$newfamily{'sable'}=\@accarray;
	$newfamily{'sablefile'}="$sabledir$id_RES";
	#my $acc=calcacc(\%hit);
	#my $accscore=acccal($acc);
	#$family{'aacfile'}=$seqdir.$id.".acc";
	#$family{'solacc'}=$accscore;
	#$family{'acc'}=$acc;
	}
	push (@outarray, \%newfamily);
	} 
	}
return @outarray;
}

sub temphithash
{
my $id=$_[0];

$id=~m/(.{4,9})\:(.)\:(\d{1,4})\-(\d{1,4})/;
my $pdbid=$1;
my $chain=$2;
my $motifstart=$3;
my $motifmid=$3+2;
my $motifend=$4;
my @res=($motifstart,$motifmid,$motifend);
my %hit=('file',$id,'chain',$chain,'pdbid',$pdbid,'resnum',\@res);
return %hit;
}

sub getfamilydata
{
my $startfile=$_[0];
my $Newquery=$_[1];
my $dir=$_[2];
my $blastdir=$_[3];
my $seqdir=$_[4];
my $blastdb=$_[5];
my $sabledir=$_[6];
my @families=parsexml($startfile);
my @outarray;
my $familycnt=@families;
#my $familycnt=1;

print "$familycnt families found\n";

for(my $i=0;$i<$familycnt;$i++)
	{
	my %family=%{$families[$i]};
	my @idarray=@{$family{pdblist}};
	#my $idcnt=@idarray;
	my $idcnt=1;
for (my $j=0; $j<$idcnt;$j++)
{	my $id=$idarray[$j];
	my %hit=temphithash($id);
	my $pdb=$hit{pdbid};
my $file="$dir$pdb.pdb";
if(open(FILE ,$file))
{
print "$file found...";
}
elsif (!$Newquery)
{
my $pdb=$hit{pdbid};
pdbftp($dir,$pdb);
}

	print "\n\nfamily $i+1 of $familycnt, $id used as representative\n";

	
	#print "running phiblast...\n";
	#phiblast(\%hit,,$Newquery,$dir,$blastdir,$seqdir);
	
	print "running psipred...\n";
	psipred	(\%hit,$Newquery,$dir,$blastdir,$seqdir);
	alignsecondary($id);

	print "calculating solvent accesability...\n";
	runsable(\%hit,$Newquery,$dir,$seqdir,$sabledir);
	}
	} 

}

sub getfamilydata2
{
my $startfile=$_[0];
my $Newquery=$_[1];
my $dir=$_[2];
my $blastdir=$_[3];
my $seqdir=$_[4];
my $blastdb=$_[5];
my $sabledir=$_[6];
my @families=parsexml($startfile);

my @outarray;
my $familycnt=@families;
#my $familycnt=1;
print "START:$startfile\nNEW:$Newquery\nDir:$dir\nBlastdir:$blastdir\nBlastDB:$blastdb\nSabledir:$sabledir\n";
print "$familycnt families found\n";

for(my $i=0;$i<$familycnt;$i++)
	{
	my %family=%{$families[$i]};
	my @idarray=@{$family{pdblist}};
	print "IDARRAY:@idarray\n";
	#my $idcnt=@idarray;
	my $idcnt=1;
for (my $j=0; $j<$idcnt;$j++)
{	my $id=$idarray[$j];
	my %hit=temphithash($id);
	my $pdb=$hit{pdbid};
print "PDB: $pdb\n";


my $seqfile="$seqdir$pdb.fas";
print "FILE:$seqfile\n";

	print "\n\nfamily $i+1 of $familycnt, $id used as representative\n";

	
	#print "running phiblast...\n";
	#phiblast(\%hit,,$Newquery,$dir,$blastdir,$seqdir);
	
	#print "running psipred...\n";
	#psipred	(\%hit,$Newquery,$dir,$blastdir,$seqdir,$seqfile);
	#alignsecondary($id);

	print "calculating solvent accesability...\n";
	runsable(\%hit,$Newquery,$dir,$seqdir,$sabledir,$seqfile);
	}
	} 

}

sub getfamilydatap
{
my @families=@{$_[0]};
my $Newquery=$_[1];
my @outarray;
my $familycnt=@families;
#my $familycnt=1;

print "$familycnt families found\n";

for(my $i=0;$i<$familycnt;$i++)
	{
	my %family=%{$families[$i]};
	my @idarray=@{$family{pdblist}};
	#my $idcnt=@idarray;
	my $idcnt=1;
for (my $j=0; $j<$idcnt;$j++)
{	my $id=$idarray[$j];
	my %hit=temphithash($id);
	my $pdb=$hit{pdbid};
my $file="/home/dwoodhead/wrkdir/ftpdir/$pdb.pdb";
if(open(FILE ,$file))
{
print "$file found...";
}
elsif (!$Newquery)
{
my $pdb=$hit{pdbid};
print "pdb not found or incomplete\n: getting:$pdb\n";
	my $localfile="/home/dwoodhead/wrkdir/ftpdir/pdb$pdb.ent.gz";
	system("wget -O $localfile ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb$pdb.ent.gz");
 	if (open(FILE,"<$localfile"))
		{
		system("uncompress /home/dwoodhead/wrkdir/ftpdir/pdb$pdb.ent.gz");
		system("cp /home/dwoodhead/wrkdir/ftpdir/pdb$pdb.ent /home/dwoodhead/wrkdir/ftpdir/$pdb.pdb") ;
		system("rm /home/dwoodhead/wrkdir/ftpdir/pdb$pdb.ent");
		}
	else 
		{
		print "file not downloaded\n";
		open (BINDLOG, ">$resultsdir.$runident.bindinglog");
		print BINDLOG "file $pdb not found in pdb database $filename not made.\n"; 
		} 
}

	print "\n\nfamily $i+1 of $familycnt, $id used as representative\n";

	
	print "running phiblast...\n";
	phiblast(\%hit,,$Newquery);
	
	print "running psipred...\n";
	psipred	(\%hit,$Newquery);
	alignsecondary($id);

	print "calculating solvent accesability...\n";
	runsable(\%hit,$Newquery);
	}
	} 

}
sub getgendata
{
my @hits=@{$_[0]};
my $Newquery=$_[1];
my $seqdir=$_[2];
my $blastdir=$_[3];
my $sabledir=$_[4];
my $phiblast=$_[5];
my $psipred=$_[6];
my $sable=$_[7];
my @outarray;
my $familycnt=@hits;
#my $familycnt=1;

print "$familycnt families found\n";

for(my $i=0;$i<$familycnt;$i++)
	{
	my %hit=%{$hits[$i]};
	#my $idcnt=@idarray;
	my $idcnt=1;
	my $id=$idarray[$i];
	my $pdb=$hit{pdbid};
        print "SEQDIR:$seqdir\nID:$pdb\n";
	$seqfile="$seqdir$pdb.fas";
	print "\n\nfamily $i+1 of $familycnt, $id used as representative\n";

	if ($phiblast==1)
	{
	print "running phiblast...\non $seqfile\n";
	phiblast(\%hit,"",$Newquery,$ftpdir,$blastdir,$seqdir,$blastdb,$seqfile);
	}
	if ($psipred==1)
	{
	print "running psipred...\non $seqfile\n";
	psipred	(\%hit,"",$blastdir,$seqdir,$seqfile);
	alignsecondary($id,$blastdir);
	}
	if ($sable==1)
	{
	print "calculating solvent accesability...\non $seqfile\n";
	runsable(\%hit,$Newquery,$seqdir,$sabledir,$seqfile);
	}
	} 

}

sub makegenseqfiles
{
my @hitlist=@{$_[0]};
my $seqdir=$_[1];
my $fastafile=$_[2];


$fastafile=~m/(.*?)\..*/;
$fastaout=$1."_scan.fas";

my $hitcnt=@hitlist;
print"making sequence files\n $hitcnt hits found\n";
#if(open(FILE, "$fastaout"))
#	{}
#	else
#	{
	print "making file $fastaout\n";
	open (FASTAOUT,">$fastaout");
for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my $id=$hit{pdbid};
	my  $sequence=getgenseq(\%hit,"","",$fastafile);
	print FASTAOUT ">$id\n$sequence";
	}
	
#	}

print"making sequence files\n $hitcnt hits found\n";

for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my $id=$hit{pdbid};
	
	if(open(FILE, "$seqdir$id.fas"))
	{}
	else
	{
	
	print "making sequence file $seqdir$id.fas\n";
	#print "SEQ:$sequence\n";
	my  $sequence=getgenseq(\%hit,"","",$fastafile);
	open (SEQFILE,">$seqdir$id.fas");
	print SEQFILE ">$id\n$sequence"; 
	}
	}

}

sub pickpositive
{
my $familyfile=$_[0];
my @families=parsexml($familyfile);
my @outarray;

my $familycnt=@families;
print "$familycnt families in array\n";

for (my $i=0;$i<$familycnt;$i++)
{
my %familydata=%{$families[$i]};
my $famid=$familydata{scop};
if($famid=~m/^(\w\.\d{1,4}\.\d{1,4})\.\d{1,4}/)
{
my $superfamilyid=$1;
$familydata{scop}=$superfamilyid;
}
else{}
$bvalue=$familydata{bvalue};
#print "$famid binding value:$bvalue\n";
if ($bvalue eq "B+")
	{
#	print "adding to out array\n";
	push(@outarray,\%familydata);
	}
else
	{
#	print "family discarded\n";
	}
}
return @outarray;

}

1;
