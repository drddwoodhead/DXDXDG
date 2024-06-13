sub testmetal #tests if metal atom is present associated with any of the residues matched
#arguments are the 4 character pdb id, the location of the pdb file and the residues required.
{
## arguments
my $pdb=$_[0];
$pdb=lc($pdb);
my $dir=$_[1];
my $chain=$_[2];
my @met=@{$_[3]};
my @res=@{$_[4]};
my $returnflag=$_[5];
my$totres=@res;
my$threshold=4;
my @candidatemet;
#print "\n\npdbid $pdb\n";
#print "directory $dir\n";
#print "@_\n";
#print " chain:$chain\n";
#print ("number of residues:$totres\n");
$openpdb ="$dir$pdb.pdb";
#print "locatoin of pdb :$openpdb\n";

$filename =$pdb.":".$chain.":".$res[0]."\-".$res[$totres-1]."\.pdb";
#print "\n\ntesting $pdb\n";
#print ("file created: $filename\n");
my@metxyz;
my@Oxyz=getOligands($openpdb,\@res, $chain);
if (!$metxyz[0])
	{
	@metxyz=gethetatoms($openpdb);
	}
my $Otot=@Oxyz;
my $mettot=@metxyz;
#print "$Otot oxygen atoms associated with residues and $mettot hetatoms found in file\n";
  
for  (my$j=0; $j<$mettot; $j++)
	{ 
my @closeres;
	my %every=%{$metxyz[$j]};
	#print "metal $every{met}\n";
	for  (my$i=0; $i<$Otot; $i++)
  		{
	
  		my %each=%{$Oxyz[$i]};
		#print "residue $each{res}\n";	

		my $distance=pythag3d(\%each,\%every);
		#print "distance is : $distance\n";
		if ($distance<=$threshold)
			{
			#print "metal ion found near binding residue: $every{met}\n";
			my $flag=0;
			#print "residues found near metal so far @closeres\n";
			#print "$distance\n";
			for $resfound(@closeres)
				{
				if($each{res} eq $resfound)
					{$flag=1;}
				else
					{}
				}
			if (!$flag)
				{
				#print "$every{met} found near res: $each{resid}:$each{res}\n";
				push(@closeres,$each{res});
				}	
	
			}
		else 
			{}
		}
my $closeres=@closeres;
undef(@readpdb);
#print "TOTres:$totres==$closeres:closeres\n";
if ($closeres==$totres)
{
#print "closeres:$closeres=$totres:totres\n";
#print"candidtatemet line : $every{line}\n";
push (@candidatemet, \%every);
}
else{}
}



#@candidatemet= sort {$$a{ave} <=> $$b{ave}} @candidatemet;

my %choosen=%{$candidatemet[0]};
if (%choosen)
{
my $count=@candidatemet;
#print"found $count metals within 3A of all 3 residues\n ";
#print "metxyz=$choosen{x},$choosen{y},$choosen{z},  $choosen{met}\n";
return @candidatemet;
#$return=$choosen{met}.",".$choosen{ave};
#return $return;
}
else
{
 return 0;
}

}

sub findligands #tests if metal atom is present associated with any of the residues matched
#arguments are the 4 character pdb id, the location of the pdb file and the residues required.
{
## arguments
my $pdb=$_[0];
$pdb=lc($pdb);
my $dir=$_[1];
my $chain=$_[2];
my @metxyz=@{$_[3]};
my @res=@{$_[4]};
my$totres=@res;
my$threshold=4;
my @candidateligands;
my @outligands;
my $sitedir="/home/dwoodhead/wrkdir/sitepdb/";
my $openpdb ="$dir$pdb.pdb";
#print "locatoin of pdb :$openpdb\n";

$opensitepdb=$sitedir.$pdb.":".$chain.":".$res[0]."\-".$res[$totres-1]."\.pdb";
#print "\n\ntesting $pdb\n";
#print ("file created: $filename\n");
my@nullres;
my@Oxyz=getOligands($openpdb);
my @metxyz=gethetatoms($opensitepdb);

  my $Otot=@Oxyz;
  my $mettot=@metxyz;
#print "$Otot oxygen atoms associated with residues and $mettot hetatoms found in file\n";

for  (my$j=0; $j<$mettot; $j++)
	{
my @closeres;
	my %every=%{$metxyz[$j]};
	#print "metal $every{met}\n";
	for  (my$i=0; $i<$Otot; $i++)
  		{
		my %each=%{$Oxyz[$i]};
		#print "residue $each{res}\n";	
		my $distance=pythag3d(\%each,\%every);
		#print "distance is : $distance\n";
		if ($distance<=$threshold)
			{
			#print "LINE: $each{line}\n";
			my %push=%each;
			push (@candidateligands, \%push);	
			}
		else 
			{}
		}
	}
my $count=@candidateligands;
#print"found $count O within $threshold A of the metal\nremoving residues associated with motif...\n";
for (my $i=0;$i<$count;$i++)
	{
	#print "opening new potential Ligand $i of $count\n";
	my %ligand=%{$candidateligands[$i]};
	#print "\nLIGANDLINE:$ligand{line}";
	my $resflag=0;
	for my $testres(@res)
		{
		#print "res:$ligand{res}:$testres:\n";
		if ($ligand{res} == $testres)
			{
			#print "res:$ligand{res}\n";
			$resflag=1;
			}
		}
	if ($resflag==0)# && ($ligand{resid} eq "ASN" ||$ligand{resid} eq "ASP"))
		{
		#print "potential new ligand:".$ligand{line};
		push (@outligands,\%ligand);
		}
	}
my $count=@outligands;
#print"found $count O within $threshold A of the metal\n ";
return @outligands;
#$return=$choosen{met}.",".$choosen{ave};
#return $return;
}

sub addliganddata #tests if metal atom is present associated with any of the residues matched
#arguments are the 4 character pdb id, the location of the pdb file and the residues required.
{
#print "searching for bound metals\n";
my @hitlist=@{$_[0]};
my $hitlist=@hitlist;
my @outlist;$openpdb ="$dir$pdb.pdb";
for  ($i=0; $i<$hitlist; $i++)
	{
	  my%each=%{$hitlist[$i]};
	
	
	my $pdbid=$each{pdbid};
	#print "pdbid:$pdbid";
        my $chain=$each{chain};
	my $pdbfileid=$each{file};
	my $dir="/home/dwoodhead/wrkdir/ftpdir/";
	#print ":$chain:@res\n";
	my $opensitepdb ="$dir$pdbfileid.pdb";
	my @res=@{$each{resnum}};
	
	my @ligand=findligands($pdbid,$dir,$chain,\@met,\@res);
	if (@ligand)
	 	{
	 	$each{ligand}=\@ligand;
		push(@outlist,\%each);
		}
	else
		{
		}
	
	}
undef(@hitlist);
return @outlist;
}

sub gethetatoms
{
my $openpdb=$_[0];
#print "opening file $openpdb to find HETATOMS\n";
my @metxyz;

open(PDBFILE, $openpdb);
 @readpdb = <PDBFILE>;
 close PDBFILE;
READLINE1:for my $line(@readpdb)
  { 
 #if( $line =~ /^HETATM.{6}\s{0,2}CA|\s{0,2}MG|\s{0,1}UNK|\s{0,2}OH/)
 if( $line =~ /^HETATM.*/)
        { #reads in coordinate of metal ions
        my $metion=substr($line,17,3); #print "metal ion is: $metion\n";
	$metion=~s/\s//;
	#my $line=~m/^HETATM.{10}\s{0,4}(\w{0,4})/; $metion=$1;
	#print "metal ion is: $metion\n";
        my $metx= substr($line,30,8); #print "mx is: $metx\n"; ##### check these coordinates
        my $mety= substr($line,38,8); #print "My is: $mety\n";
        my $metz= substr($line,46,8); #print "Mz is: $metz\n\n";
	my $thisline=$line;
	#print "get hetatoms line:$thisline\n";
        my%metcoordinates=('met',$metion,'x',$metx,'y',$mety,'z',$metz,'line',$thisline);
        #print"$metion\n";
	#print "metxyz=$choosen{x},$choosen{y},$choosen{z},  $choosen{met}\n";
        push(@metxyz, \%metcoordinates );
	} 
  elsif($line =~ /^\s{0,4}ENDMDL/)
	{
	#print PDBOUT "END";
	last READLINE1
	}
     elsif($line =~ /^\s{0,4}END/)
	{
	#print PDBOUT "END";
	last READLINE1
	}
      else{}
  }
return @metxyz;
}

sub getOligands
{

my $openpdb=$_[0];
#print "opening file $openpdb to find O\n";
my @res=@{$_[1]};
my $chain=$_[2];
my @Oxyz;
my $flag;
#print "RES:@res\n";
open(PDBFILE, $openpdb);
 @readpdb = <PDBFILE>;
 close PDBFILE;
if
($res[0]==0)
{
$res[0]="runonce";
#print"no residues supplied\n";
$flag=1;
}
 for my $each (@res)
 {
 #print "\nLooking at residue :$each:\n";

READLINE2:for my $line(@readpdb)
  {
	my $thischain=$chain;
	if ($thischain eq "-")
		{$thischain=" ";}
     	if( $line =~ /^ATOM.{9}[ON][GDE].{6}$thischain\s{0,3}$each/)
        { #reads in coordinate of oxygens in D or E
        my $Ox= substr($line,30,8); 
        my $Oy= substr($line,38,8); #print "Oy is: $Oy\n";
        my $Oz= substr($line,46,8); #print "Oz is: $Oz\n\n";
	my $resid= substr($line,17,3);
	
	$resid=~s/\s//;
        my %Ocoordinates=('x',$Ox,'y',$Oy,'z',$Oz,'res',$each,'resid',$resid);
	#print $line;
	#print "res:$res Oxyz=$Ocoordinates{x},$Ocoordinates{y},$Ocoordinates{z}\n";
	push(@Oxyz, \%Ocoordinates);
        }
	elsif($line =~ /^ATOM.{9}[ON][GDE]/&& $flag)
        { #reads in coordinate of oxygens in D or E
        my $Ox= substr($line,30,8); 
        my $Oy= substr($line,38,8); #print "Oy is: $Oy\n";
        my $Oz= substr($line,46,8); #print "Oz is: $Oz\n\n";
	my $resid= substr($line,17,3);
	my $resnum=substr($line,22,4);
	$resnum=~s/\s//g;
        my%Ocoordinates=('x',$Ox,'y',$Oy,'z',$Oz,'res',$resnum,'resid',$resid, 'line', $line);
	#print "Oxyz=$Ocoordinates{x},$Ocoordinates{y},$Ocoordinates{z}\n$line\n";
	push(@Oxyz, \%Ocoordinates);
        }
  	elsif($line =~ /^\s{0,4}ENDMDL/)
	{
	#print PDBOUT "END";
	last READLINE2
	}
    	 elsif($line =~ /^\s{0,4}END/)
	{
	#print PDBOUT "END";
	last READLINE2
	}
     else{}
  }
 }
return @Oxyz
}

sub pythag3d #takes 2 hash tables with 3d coordinates and returns their 3d seperation.
{
%p1=%{$_[0]};
%p2=%{$_[1]};

$l=$p1{'x'}-$p2{'x'};
$h=$p1{'y'}-$p2{'y'};
$d=$p1{'z'}-$p2{'z'};

$sqdistance=($l*$l)+($h*$h)+($d*$d);
#print "$p1{'x'},$p1{'y'},$p1{'z'}$p2{'x'}$p2{'y'}$p2{'z'}\n";

$distance= sqrt($sqdistance);

return $distance;

}

sub addmetalbinding
{
#print "searching for bound metals\n";
my @hitlist=@{$_[0]};
my $ftpdir=$_[1];
my $sitedir=$_[2];
my $hitlist=@hitlist;
my @outlist;
#print "sitedir***:$sitedir\n";
for  ($i=0; $i<$hitlist; $i++)
	{
	my%each=%{$hitlist[$i]};
	my $pdbid=$each{pdbid};
#	print "ID :pdbid:$pdbid\n";
	my @res=@{$each{resnum}};
        my $chain=$each{chain};
	my $pdbfileid=$each{file};
	my $restot=@res;
	my $lastres=$res[$restot-1];
	
	
#	print ":$chain:@res\n";
	my @met;
	my $openpdb="$ftpdir$pdbid.pdb";
	my $opensitepdb=$sitedir.$pdbid."_".$chain."_".$res[0]."-".$lastres.".pdb";
	if (!open(PDBFILE,"<$openpdb"))
		{ print "pdb not found or incomplete at $openpdb\n";
		}
	my @metal=testmetal($pdbid,$ftpdir,$chain,\@met,\@res);
	my $metcnt=@metal;
	if (%{$metal[0]}==0){$metcnt=0;}
	#print "$metcnt metals found\n";
	#print "METAL:".%{$metal[0]}."\n";
	if ($metcnt != 0)
		{
		$each{metal}=\@metal;
		my $metcount=@metal;
	 
	#print "opening file $opensitepdb\n";
	open(READSITEPDB,"<$opensitepdb") || die("Cannot Open File"); 
	@readsitepdb=<READSITEPDB>;
	close READSITEPDB;
	open(PDBOUT,">$opensitepdb") || die("Cannot Open File"); 
	
	for my $line(@readsitepdb)
 		{
     		if( $line =~ /^\s{0,4}ATOM/)
        		{
        		print PDBOUT "$line";
        		}
     		else {}
  		}
	for (my $i=0;$i<$metcount;$i++)
		{
		my%met=%{$metal[$i]};
		my $metline=$met{line};
	#	print "METLINE:$metline\n";
		print PDBOUT $metline;
		}
	print PDBOUT "END";
		
	close PDBOUT;
	}
	 else{}
	 push(@outlist,\%each);
	}
undef(@hitlist);
	return @outlist;
}

sub findbound
{
my @hitlist=@_;
my @outarray;
my$hitlistcnt=@hitlist;

for (my $i=0;$hitlistcnt>$i;$i++)
{
my %hit=%{$hitlist[$i]};
#print "metal bound:$hit{metal}\n";


if($hit{metal})
{
#print "metal bound :true\n";
push(@outarray,\%hit);
}

}
my $count = @outarray;
#print "$count hits left...\n";
return @outarray;
}

sub sortbybinding
{
my @bindinglist;

my @hitlist=@{$_[0]};
$hitlist=@hitlist;

for  (my $x=0 ; $x<$hitlist; $x++)
	{
	my %each=%{$hitlist[$x]};
	my $pdbid=$each{file};
	my @metal=@{$each{metal}};
	my $metcnt=@metal;
	for (my $j=0; $j<$metcnt;$j++)
		{
		my %metal=%{$metal[$j]};
		my $metal=$metal{met};
		my $bindingflag=0;
		my %selectedbinding=0;
		my $arrayref=0;

		my$bindinglist=@bindinglist;	
	
		for  (my $j=0 ; $j<$bindinglist; $j++)
			{
			my%bindingentry=%{$bindinglist[$j]};
			my$recordedmetal=$bindingentry{metal};
			#print"i=$x j=$j comparison $metal == $recordedmetal\n";
			my $metcnt=@metal;
	
			if($metal eq $recordedmetal)
		 		{
          			#print "TRUE setting flag to 1\n";
	  			$bindingflag=1;
	  			%selectedbinding=%bindingentry;
	  			$arrayref=$j;	
		 		}
			else 	
				{
				#print "flag is $bindingflag\n";
				}
			
			}
		if ($bindingflag)
			{
			my@pdblist=@{$selectedbinding{pdblist}};
			push(@pdblist,$pdbid);
			$selectedbinding{pdblist}=\@pdblist;
			$bindinglist[$arrayref]=\%selectedbinding;
			#print "binding true adding $pdbid to binding $selectedbinding{metal}\n";
			}
		else
			{
			my%newbinding;
			$newbinding{metal}=$metal;
			my@pdblist;
			push(@pdblist,$pdbid);
			$newbinding{pdblist}=\@pdblist;
			#$newbinding{binding}=;#placeholder for the written bindingification.
			#print "new binding found, creating binding $metal, with pdb $pdbid\n"; 

			push(@bindinglist,\%newbinding)
			}
		}
	}	

return @bindinglist;
}

sub isbound
{
 my @hitlist=@{$_[0]};
 my $ID=$_[1];
my $hitlistcnt=@hitlist;
for (my $i=0; $i<$hitlistcnt; $i++)
{
my %hit=%{$hitlist[$i]};
my $hitid=$hit{file};
if ($ID eq $hitid)
{
return $hit{bound};
}
else{}
}
}

sub pdbftp
{
my $pdbid=$_[1];
my $ftpdir=$_[0];
my $open=$ftpdir.$pdbid.".pdb";
if (open(PDB,"<$open"))
	{
	print"$ftpdir$pdbid found\n";
	}
else
{
print "getting:$pdbid\n";
		my $localfile=$ftpdir."pdb$pdbid.ent.gz";
		my $localfile2=$ftpdir."pdb$pdbid.ent";
		system("wget -q -O \"$localfile\" ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb$pdbid.ent.gz");

 		if (open(PDBGZ,"$localfile"))
			{
			system("uncompress \"$ftpdir"."pdb$pdbid.ent.gz\"");
			}
			if (open(PDBGZ2,"$localfile2"))
			{
			system("cp \"$ftpdir"."pdb$pdbid.ent\" \"$ftpdir$pdbid.pdb\"") ;
			system("rm \"$ftpdir"."pdb$pdbid.ent\"");
			}
		else 
			{print "file not downloaded\n";}
}
}
1;