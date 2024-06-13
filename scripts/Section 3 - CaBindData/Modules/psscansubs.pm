use clanssubs; 

sub parsepsout
{
my $file=$_[0];
my $gen=$_[1];
my @hitlist;
open(PSFILE,"<$file");
my @psfile=<PSFILE>;
close PSFILE;
my $pscnt=@psfile;

print "openfile: $file...\n$pscnt lines read...\n";

for(my $i=0;$i<$pscnt;$i++)
	{
	my $line=$psfile[$i];
	if ($line=~m/^\>(.{4})([\w\s])\s/)
		{
		$id=$1;
		$id=lc($id);
		$chain=$2;
		#print "ID:$id, Chain:$chain\n";
		}
	elsif($line=~m/^\>(.*?)\s/)
		{
		$id=$1;
		$id=lc($id);
		$chain="-";
		#print "ID:$id.\n";
		}
	if ($line=~m/^\s{3,7}(\d{1,5})\s\-\s(\d{1,5})\s\s(\w*)/)
		{
		
		my $res1=$1;
		my $res2=$res1+2;
		my $res3=$2;
		my $motif=$3;
		$motif=uc($motif);
		if (!$gen==1)
		{
		my $seqstart=getstartresidue($openpdb,$chain);
		if ($seqstart==0){$res1--;$res2--;$res3--;} #hack as ps scan seems to incorrectly report residue numbers if first residue in chain is listed as 0
		}
		my $fileid=$id.":".$chain.":".$res1."-".$res3;
		my @res=($res1, $res2, $res3);
		print "$id,$fileid,$chain,$res1,$res2,$res3,$motif\n";
		my %hit=('pdbid',$id,'file',$fileid,'chain',$chain,'resnum',\@res,'motif',$motif);
		push (@hitlist, \%hit);
		}
	}
return @hitlist;
}

sub downloadpdb
{

my $file="/home/dwoodhead/wrkdir/ftpdir/$id.pdb";
if(open(FILE ,$file))
{
print "$file found...";
}
else
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
		{print "file not downloaded\n";
		open (BINDLOG, ">$resultsdir.$runident.bindinglog");
		pint BINDLOG "file $pdb not found in pdb database $filename not made.\n"; 
		} 
}

}

sub runpsscan
{
my $seqfile=$_[0];
my $psscanfile=$_[1];
my $pattern=$_[2];

open (FILE ,"<$seqfile");
my @file=<FILE>;
close FILE;
my @file2;
open (FILE, ">$seqfile");
for $each(@file)
{
if ($each=~m/^>gi\|(.*?)\|.*/)
{
#print "$1\n";
$each=">".$1."\n";
print FILE $each;
}
elsif ($each=~m/^(>.*?)\|(.*)\|.*/)
{
#print "$1\n";
$each=">".$2."\n";
print FILE $each;
}
else 
{
$each=~s/ /_/g;
print FILE $each;
}
}
system(`perl /home/dwoodhead/ps_scan/ps_scan.pl -p $pattern "$seqfile">"$psscanfile"`);
#system(`perl /home/dwoodhead/ps_scan/ps_scan.pl -p D-x-[DNS]-x[DNS] "$seqfile">"$psscanfile"`);
}


1;