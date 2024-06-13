sub genomedata
{
my @hitlist=@{$_[0]};
my $genomefile=$_[1];

my $hitcnt=@hitlist;
my @outlist;
#print "genomedata: $hitcnt Families found\n";
#print "genomefile : $genomefile\n";
for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my $id=$hit{repid};
		my $gendata=getgendata($genomefile,$id);
	print "ID:$id:$gendata\n";
	$hit{gendata}=$gendata;
	push(@outlist,\%hit);
	}
return @outlist;
}

sub getgendata
{
my $fastafile=$_[0];
my $seqid=$_[1];

$seqid=~m/(\d{8,9})\:\-\:.*/;
#print "1:$1,2:$2,3:$3,4:$4,5:$5,6:$6,7:$7\n";
my $id=$1;
#print "ID:$id\n";

open(FASTAFILE,"<$fastafile");
@fastafile=<FASTAFILE>;
close FASTAFILE;

$linecnt=@fastafile;
my $gendata;

for(my $i=0;$i<$linecnt;$i++)
{

my $line=$fastafile[$i];
if ($line=~m/>gi\|$id\|(.*)/)
	{
	$gendata=$1;
	#print "$seqid:$gendata\n";
	}

else
	{
	#print "false\n"
	}
}
return $gendata;
}

sub getrpsdata
{
my @rpsfile=@{$_[0]};
my $id=$_[1];
my %rpshash=%{$_[2]};
my $rpsfile=$_[3];



$linecnt=@rpsfile;
#print "$rpsfile :$linecnt Lines\n";
my @rpsdata;


#print "ID:$id\n";



my $index=$rpshash{$id};

#print "ID $id: LINE $index\n";

for(my $i=$index;$i<$linecnt;$i++)
{
my $line=$rpsfile[$i];
            #gnl|CDD|109572 pfam00521, DNA_topoisoIV, DNA gyrase/topoisomeras...   757   0.0  

if ($line=~/^Query\=(.*)/&& $inquery==0)
	{
	$qid=$1;
	$qid=~m/(\w\.\d{1,4}\.\d\.\d)\:(.*?)\:(\w)\:(\d{1,4})-(\d{1,4})/;
	my @res=($4,$4+2,$5);
	my $id="$2:$3:$4-$5";
	%newquery=('id',$id,'scop',$1,);
	}
elsif($line=~m/gnl\|CDD\|\d{1,7}\s(.*?),(.*?\.)\s{1,5}(\d{1,5})\s{3}(.*?)\s/)
	{
	print "MotifID:$1,Description:$2,Score:$3,E-Val:$4\n";
	
	my %newhit=('motifid',$1,'description',$2,'score',$3,'eval',$4);
	push (@tempdata,\%newhit);
	
	}
elsif($line=~m/>gnl\|CDD\|\d{1,7}\s(.*?),(.*?\.)\s/)
	{
	$inquery=0;
	$hits=@tempdata;
	for(my$i=0;$i<$hits;$i++)
		{
		
		}
	#print "end query $hits hits found\n"; 
	}
elsif($line=~/^Query\=(.*)/&& $inquery==1)
	{
	$inquery=0;
	}
else
	{
	#print "false\n"
	
}
}
}

sub rpsdata
{
my @rpsfile=@{$_[0]};
my $id=$_[1];
my %rpshash=%{$_[2]};
my $rpsfile=$_[3];
$linecnt=@rpsfile;
print "$rpsfile :$linecnt Lines\n";
my @rpsdata;
print "ID:$id\n";
my $index=$rpshash{$id};
print "index:$index\n";
for(my $i=$index;$i<$linecnt;$i++)
{
my $line=$rpsfile[$i];
            #gnl|CDD|109572 pfam00521, DNA_topoisoIV, DNA gyrase/topoisomeras...   757   0.0  
#print "$i:$line";
if($line=~m/gnl\|CDD\|\d{1,7}\s(.*?),(.*?\.)\s{1,5}(\d{1,5})\s{3}(.*?)\s/)
	{
	print "*rpsdata:MotifID:$1,Description:$2,Score:$3,E-Val:$4\n";
	
	my %newhit=('motifid',$1,'description',$2,'score',$3,'eval',$4);
	push (@rpsdata,\%newhit);
	
	}
elsif($line=~m/>gnl.*/)
	{
	$inquery=0;
	$hits=@rpsdata;
#	print "end query $hits hits found\n"; 
	return @rpsdata;
	}

else
	{
	#print "false\n"
	}
}
}
sub rpsblastdata
{
my @hitlist=@{$_[0]};
my $rpsfile=$_[1];

my $hitcnt=@hitlist;
#my $hitcnt=10;
my @outlist;
print "rpsdata: $hitcnt Families found\n";
print "rpsfile : $rpsfile\n";

open(RPSFILE,"<$rpsfile");
@rpsfile=<RPSFILE>;
close RPSFILE;

my %rpshash=rpshash(\@rpsfile);
my $lastid;
my $lastrpsdata;

for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my $id=$hit{repid};
	$id=~m/(\d{8,9})\:\-\:.*/;
	$id=$1;
	if($id==$lastid)
	{
	$hit{rpsdata}=\@lastrpsdata;
	push(@outlist,\%hit);
	}
	else
	{
	print "rpsblastdata $i of $hitcnt ID:$id\n";
	
	my @rpsdata=rpsdata(\@rpsfile,$id,\%rpshash,$rpsfile);
	$hit{rpsdata}=\@rpsdata;
	push(@outlist,\%hit);
	$lastid=$id;
	@lastrps=@rpsdata;
	}
	}
return @outlist;
}

sub rpsblastdata2
{
my @hitlist=@{$_[0]};
my $rpsfile=$_[1];

my $hitcnt=@hitlist;
my @outlist;
#print "rpsdata: $hitcnt Families found\n";
#print "rpsfile : $rpsfile\n";

open(RPSFILE,"<$rpsfile");
@rpsfile=<RPSFILE>;
close RPSFILE;

my %rpshash=rpshash(\@rpsfile);
my $lastid;
my $lastrpsdata;

for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my @idlist=@{$hit{pdblist}};
	my $idcnt=@idlist;
#	print "$idcnt pdbs in family\n";
	for (my $j=0;$j<$idcnt;$j++)
	{
	
	my $id=$idlist[$j];
	#print "$id,$lastid\n";
	if($id eq $lastid)
	{
	#print"$id eq $lastid\n";
	$hit{rpsdata}=\@lastrpsdata;
	push(@outlist,\%hit);
	}
	else
	{
	print "rpsblastdata2 $i of $hitcnt, $j of $idcnt, ID:$id\n";
	my @rpsdata=rpsdata(\@rpsfile,$id,\%rpshash,$rpsfile);
	my %newhit=%hit;
	$newhit{rpsdata}=\@rpsdata;
	$newhit{id}=$id;
	push(@outlist,\%newhit);
	$lastid=$id;
	@lastrps=@rpsdata;
	}
	}
	}
return @outlist;
}

sub rpshash
{
my @search=@{$_[0]};
my $cnt=@search;
my %outhash;

print "making hash table of rps file.";


#print "RPSFILE: @search\n";

for (my $i=0; $i<$cnt; $i++)
{
my $line=$search[$i];

if ($line=~m/^Query=\s(.*)/)
{
my $id=$1;
$outhash{$id}=$i;
#print "ID:$id Line:$i\n";
}
else
{}
}
return %outhash;
}

sub parsecabindresults
{
my $openfile=$_[0];
my @currentlist=@{$_[1]};
my $flag=$_[2];


open(OPENFILE,"<$openfile");
my @openfile=<OPENFILE>;
close OPENFILE;

my $linecnt=@openfile;
print "open file:$openfile\n line count:$linecnt\n";
for (my $i=0;$i<$linecnt;$i++)
	{
	my $line=$openfile[$i];
	if ($line=~m/Results:(.*)/)
		{
		$resultid=$1;
		#print "LINE:$line";
		#print "Resultid:$resultid\n";
		}
	elsif($line=~m/Positive\t(.*?)(\t.*)/||$line=~m/True\sPositive\t(.*?)(\t.*)/||$line=~m/False\sPositive\t(.*?)(\t.*)/)
		{
		#print "LINE:$line";
		my $positivecnt=$1;
		#print "positive count:$positivecnt\n";
		my $idlist=$2;
		my @idlist=split(/\t/,$idlist);
		my $idcnt=@idlist;
		#print "IDcount:$idcnt\n";
		for(my $i=0; $i<$idcnt; $i++)
			{
			my $id=$idlist[$i];
			$id=~s/\s//g;
			my $listcnt=@currentlist;
			my $idfound=0;
			#print "LISTCOUNT:$listcnt\n";
			for (my$j=0;$j<$listcnt;$j++)
				{
				my %current=%{$currentlist[$j]};
				my $curid=$current{repid};
				if ($id eq $curid)
					{
					#print "$curid=$curid\n";
					$idfound=1;
					push(@{$current{positive}} ,$resultid);
					$currentlist[$j]=\%current;
					}
				else
					{
					}
				}
			if(!$idfound)	
				{
				my @positive=($resultid);
				my %new=('repid',$id,'positive',\@positive);
				if($line=~m/True\sPositive\t.*/)
				{
				print "known value: B+\n";
				$new{knownvalue}='B+';
				}
				elsif($line=~m/False\sPositive\t.*/)
				{
				print "known value: B-\n";
				$new{knownvalue}='B-';
				}
				push(@currentlist,\%new);	
				}			
			}
		}
	elsif(($line=~m/Negative\t(.*?)(\t.*)/ && $flag eq "N") ||($line=~m/True\sNegative\t(.*?)(\t.*)/ && $flag eq "N") ||($line=~m/False\sNegative\t(.*?)(\t.*)/ && $flag eq "N"))
		{
		my $negativecnt=$1;
		my $idlist=$2;
		my @idlist=split(/\t/,$idlist);
		my $idcnt=@idlist;
		print "NEG IDcount:$idcnt\n";
		for(my $i=0; $i<$idcnt; $i++)
			{
			my $id=$idlist[$i];
			$id=~s/\s//g;
			my $listcnt=@currentlist;
			my $idfound=0;
			#print "LISTCOUNT:$listcnt\n";
			for (my$j=0;$j<$listcnt;$j++)
				{
				my %current=%{$currentlist[$j]};
				my $curid=$current{repid};
				if ($id eq $curid)
					{
					$idfound=1;
					}
				else
					{
					}
				}
			if(!$idfound)	
				{
				my @positive;
				my %new=('repid',$id,'positive',\@positive);
				if($line=~m/True\sNegative\t.*/)
				{
				print "known value: B-\n";
				$new{knownvalue}='B-';
				}
			elsif($line=~m/False\sNegative\t.*/)
				{
				print "known value: B+\n";
				$new{knownvalue}='B+';
				}
				push(@currentlist,\%new);	
				}			
			}
		}
	elsif($line=~m/Errors\t(\d{0,5})(\t.*)/)
		{
		my $errorcnt=$1;
		}
	}
return @currentlist;
}
sub printcollateddata
{
my @list=@{$_[0]};
my @typeid;
my $outfile=$_[1];
my @typeid=@{$_[2]};

open(OUTFILE,">$outfile");

print "Creating file: $outfile\n";
#print OUTFILE "$outfile\n"; 
my$typeidcnt=@typeid;
my $listcnt=@list;
print "LISTCNT:$listcnt\n";

for(my$i=1;$i<$listcnt;$i++)
	{
	my %hit=%{$list[$i]};
	my $id=$hit{repid};
	my $knwn=$hit{knownvalue};
	print OUTFILE "$id";
	if($knwn)
	{
	print OUTFILE "	$knwn";	
	}
	my @positive=@{$hit{positive}};
	my $poscnt=@positive;
	#print "$id:POSCNT:$poscnt\n";
	my $count=0;
	for($j=0;$j<$typeidcnt;$j++)
		{
		my $found=0;
		my $typeid=$typeid[$j];
		for($k=0;$k<$poscnt;$k++)
			{
			
			my $posid=$positive[$k];
			#print "posid:$posid,typeid:$typeid\n";
			if ($posid eq $typeid)
				{	
				$found=1;
				}
			
			}
		if ($found)
			{
			++$count;
			print OUTFILE "	$typeid";
			}
			else
			{
			print OUTFILE "	";
			}
		}
	print OUTFILE "	$count";
	my $gendata=$hit{gendata};
	print OUTFILE "	$gendata";
	print OUTFILE "\n";
	}

}

sub printXrefdata
{
my @list=@{$_[0]};
my @typeid;
my $outfile=$_[1];
my @typeid=@{$_[2]};
my %rpstable=%{$_[3]};

my @printhash=%rpstable;
print "@printhash";

open(OUTFILE,">$outfile");

print "Creating file: $outfile\n";
#print OUTFILE "$outfile\n"; 
my$typeidcnt=@typeid;
my $listcnt=@list;
print "LISTCNT:$listcnt\n";

for(my$i=1;$i<$listcnt;$i++)
	{
	my %hit=%{$list[$i]};
	my $id=$hit{repid};
	my @rpsdata=@{$hit{rpsdata}};
	my $rpsdatacnt=@rpsdata;
	my $knwn=$hit{knownvalue};
	print OUTFILE "$id";
	if($knwn)
	{
	print OUTFILE "	$knwn";	
	}
	my @positive=@{$hit{positive}};
	my $poscnt=@positive;
	#print "$id:POSCNT:$poscnt\n";
	my $count=0;
	for($j=0;$j<$typeidcnt;$j++)
		{
		my $found=0;
		my $typeid=$typeid[$j];
		for($k=0;$k<$poscnt;$k++)
			{
			
			my $posid=$positive[$k];
			if ($posid eq $typeid)
				{	
				$found=1;
				
				}
			
			}
		if ($found)
			{
			$count++;
			#print OUTFILE "	$typeid";
			}
			else
			{
			#print OUTFILE "	";
			}
			
		}
	print OUTFILE "	$count";
	print OUTFILE "\n";
	#print "rpsdatacnt:$rpsdatacnt\n";
	for(my$i=0;$i<$rpsdatacnt;$i++)
		{
		my %rps=%{$rpsdata[$i]};
		
		my $motifid=$rps{motifid};
		my @rpsref=@{$rpstable{$motifid}};
		$rpstablecnt=@rpstable;
		print OUTFILE "		$motifid\n";
#		print "ID:$motifid\n";
#		print "RPSREF: @rpsref\n";
		if (@rpsref)
		{
		#my $rpsrefcnt=@rpsref;
#		print "RPSREF: @rpsref";
		for(my$j=0;$j<$rpsrefcnt;$j++)
			{
			my %hit=%{$rpsref[$j]};
			my $id=$hit{pdbid};
			my $scop=$hit{scop};
			print OUTFILE "			$id	$scop\n";
			}
		}
		else
			{
			print OUTFILE "			no motif matches in training data\n"
			}
		}
	print OUTFILE "\n";
	}

}

sub rpsposdata
{
my @hitlist=@{$_[0]};
my $rpsfile=$_[1];


my $hitcnt=@hitlist;
my @outlist;
print "rpsdata: $hitcnt Families found\n";
print "rpsfile : $rpsfile\n";

open(RPSFILE,"<$rpsfile");
@rpsfile=<RPSFILE>;
close RPSFILE;

my %rpshash=rpshash(\@rpsfile);
my $lastid;
my $lastrpsdata;

for (my $i=0;$i<$hitcnt;$i++)
	{
	my %hit=%{$hitlist[$i]};
	my $id=$hit{id};
	
	if($id==$lastid)
	{
	$hit{rpsdata}=\@lastrpsdata;
	push(@outlist,\%hit);
	}
	else
	{
	#print "ID:$id\n";
	my @rpsdata=@{$hit{rpsdata}};
	my @rpsdata=addrpsposdata(\@rpsfile,$id,\%rpshash,\@rpsdata);
	$hit{rpsdata}=\@rpsdata;
	push(@outlist,\%hit);
	$lastid=$id;
	@lastrps=@rpsdata;
	}
	}
return @outlist;
}

sub addrpsposdata
{
my $rpsfile=$_[0];
my $id=$_[1];
my %rpshash=%{$_[2]};
my @rpslist=@{$_[3]};

my $startline=$rpshash{$id};

my $rpscnt=@rpslist;

my @outlist;

open (RPSFILE, "<$rpsfile");
my @rpsfile=<RPSFILE>;
my $filecnt=@rpsfile;

for (my $i=0;$i<$rpscnt;$i++)
	{
	my %motif=%{$rpslist[$i]};
	$motifid=$motif{motifid};
	my @positionlist;
	my $inmatch=0;	
	for (my $j=$startline; $j<$filecnt; $j++)
		{
		my $line=$rpsfile[$j];
		
		if($line=~m/^\>gnl.*/ && $inmatch)
			{
			my $inmatch=0;
			$motif{positionlist}=\@positionlist;
			last;
			}
		elsif($line=~m/^>gnl\|CDD\|\d{1,7}\s$motifid\,.*/ && !$inmatch)
			{
			$inmatch=1;
			}
		elsif($line=~m/^ Query:\s(\d{1,4})\s\w*\s(\d{1,4})/ && $inmatch)
			{
			my %pos=('qstart',$1,'qend',$2);
			push(@positionlist,\%pos);
			}
		}
	$motif{poslist}=\@positionlist;
	push(@outlist, \%motif);	

	}
return @outlist;
}

sub sortrpsdata
{
my @familydata=@{$_[0]};

my %pfamhash;

my $famcnt=@familydata;
print "sorting rps data...\n";
for (my $i=0;$i<$famcnt;$i++)
	{
	my %hit=%{$familydata[$i]};
	my @pdbidlist=@{$hit{@pdblist}};
	my $pdbcnt=@pdbidlist;
	my $scop=$hit{scop};
	my $pdbid=$hit{id};
	print "$pdbid, $scop\n";
	my @rpsdata=@{$hit{rpsdata}};
	my $rpscnt=@rpsdata;
	for (my $j=0;$j<$rpscnt;$j++)
		{
		my %rps=%{$rpsdata[$j]};
		$rps{pdbid}=$pdbid;
		$rps{scop}=$scop;
		my $motifid=$rps{motifid};
#		print "motifid:$motifid\n";
		if($pfamhash{$motifid})
			{
			my @new=@{$pfamhash{$motifid}};
			push(@new,\%rps);
			$pfamhash{$motifid}=\@new;
			}
		else
			{
			my @new=(\%rps);
			$pfamhash{$motifid}=\@new;
			}
		}	
	}
#my @printhash=%pfamhash;
#print "@printhash";
return %pfamhash;
}



1;
