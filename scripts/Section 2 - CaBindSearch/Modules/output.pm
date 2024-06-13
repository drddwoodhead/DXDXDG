use lib "/home/dwoodhead/CaBindSearch/Modules/";
use bindingsubs;

sub makelogfile
{

}

sub saveresults #makes a basic results file containg all the hits
{
 my $filename="/home/dwoodhead/wrkdir/results/"."hits\-".$seedpdb."\-".$dbfile."\-".$rmscutoff."\-".$iterations."\.smh";
 #my$filename="results:$seedpdb\_$iterations\_$dbfile";

 #print ("file created: $filename\n");
 open(RESULTS,">$filename"); # open new file

my@results=@_;
 for (my $i=0; $results[$i]; $i++)
  {
  my%each=%{$results[$i]};
  my@resnum=@{$each{resnum}};
  my@restype=@{$each{restype}};

  #print %each;
  print RESULTS "PDBID: $each{pdbid} $each{chain}\n";
 # print RESULTS "ITER: $each{iter}\n";
  %matchhash=%{$each{match}};
  print RESULTS "MATCH: %matchhash \n";
  my $restot =@{$each{resnum}};
  my $thishitinfo=$each{pdbid}.$each{resnum}[0]."\-".$each{resnum}[$restot-1];
  print RESULTS "this hit:$thishitinfo\n";
  for  (my $j=0; $results[$j]; $j++)
           	{
		my%hit=%{$results[$j]};
                my $restotc=@{$hit{resnum}};
                my $hitinfo=$hit{pdbid}.$hit{resnum}[0]."\-".$hit{resnum}[$restotc-1];

                #print RESULTS "hitinfo= $hitinfo\n";
                if ($thishitinfo eq $hitinfo)   #this hit and the hit compared is the same so rmsis 0
                 {
		 print RESULTS "RMS self ".$hitinfo.":0\n";
		 }
                if($each{$hitinfo})    #there is rms information for this hit against hit compared to
                 {
		 print RESULTS "RMS ".$hitinfo.":".$each{$hitinfo}{match}."\n";
		 }
                else              # if no rms info is found nothing to print.
                 {
		 #print RESULTS " no rms info found for $hitinfo\n"
		 }
		 #print RESULTS "---------------------------\n"
                 }
 # print RESULTS "RMS: $each{'2b59135-142'}\n";

  print RESULTS "binding: $each{metal}\n" ;
  
   for($r=0; $resnum[$r]; $r++)
    {
     print RESULTS "RES:$restype[$r] ";
     print RESULTS "$resnum[$r] : ";
     print RESULTS "$mrestype[$r]";
     print RESULTS "$mresnum[$r]\n";
    }
   print RESULTS "\n";
  }
undef(@results);
 close RESULTS;
}

sub checksitepdb
{
my $id=$_[0];
my $pdbdb=$_[1];
$id=~m/(.*?):(.):(\d{1,5})-(\d{1,5})/;
my $pdbcode=$1;
my $chain=$2;
my $startres=$3;
my $endres=$4;
my @arrayres=($3,$3+2,$4);
my $fileid=$pdbcode."_".$chain."_".$startres."-".$endres;
my $openpdb ="$pdbdb$fileid.pdb";
$sitefilecheck=0;

if (open(PDB,"<$openpdb"))
	{ 
       #if present checks pdb file is complete
	#open(PDB,"<$openpdb");
	#print "opened\n";
	my @read=<PDB>;
		for my $line(@read)
		{
			
			if ($line=~/^END/)
 				{
				#print "$line";
				#print "end of file found\n";
				$sitefilecheck=1;
				} 
			else 
				{	
				}
		}
	close PDB;
	}
else
	{ #checks pdb file is present
	#print "not opened\n"
	}


if ($sitefilecheck==1)
{
print "site pdb:$openpdb ok!\n";
}
elsif ($sitefilecheck==0)
{
print "\n$id site pdb not found or incomplete making new pdb:\n$openpdb\n";
#print "input: $pdbcode, $ftpdir, @arrayres, $chain\n";
mkpdb($pdbcode, $ftpdir, \@arrayres, $chain);
}

}

sub mkpdb #takes downloaded pdb file and creates spasm input pdb
 #arguments are the 4 character pdb id, the location of the pdb file and the residues required.
 {

 ## arguments
 my $pdb=$_[0];
 $pdb=lc($pdb);
 my $dir=$_[1];
 my $chain=$_[3];
 my $sitedir=$_[4];
my $noftp=$_[5];
#if ($chain eq "-")
#{$chain=" ";}
 my @res=@{$_[2]};
 my$totres=@res;


$openpdb ="$dir$pdb.pdb";
print "location of pdb :$openpdb\n";
my $forceftp=0;
my $filename=$sitedir.$pdb."_".$chain."_"."$res[0]\-$res[$totres-1]\.pdb";
print "making pdb file $filename\n";

if (open(FILE,"<$filename"))
{}
else
{

if (!open(FILE,"<$openpdb"))
	{ #checks pdb file is present
	$forceftp=1;
	print "$openpdb file not found\n";
	}
elsif ($noftp==0)
	{ #if present checks pdb file is complete
	open(PDB,"<$openpdb");
	@read=<PDB>;
	close PDB;	
		for my $line(@read)
		{
			$line=~ s/^ //;
			if ($line =~/^END\s.*/)
 				{
				#print "$line";
				#print "end of file found\n";
				$forceftp=0;
				} 
			else 
				{	
				$forceftp=1;
				}
		}
	}

if ($forceftp==1)
	{
	print "pdb not found or incomplete at $openpdb\n: getting:$pdb\n";
	pdbftp($dir,$pdb);
	}

open(PDBFILE, $openpdb);
 @readpdb = <PDBFILE>;
 close PDBFILE;


 #print ("file created: $filename\n");
#print "MAKING PDB : ARGUMENTS:\n";
#print "pdbid $pdb\n";
#print "directory $dir\n";
#print "@_\n";
#print " chain:$chain\n";
#print ("number of residues:$totres\n");
 open(PDBOUT,">$filename"); # open new pdb file

OUTER: for my $each (@res)
 {
  READLINE: for $line(@readpdb)
  {
     if( $line =~ /^\s{0,4}ATOM.{17}$chain\s{0,3}$each\b/||$line =~ /^\s{0,4}ATOM.{17}\s{1,4}$each\b/)
        {
        print PDBOUT "$line";
        }
     elsif( $line =~ /^\s{0,4}ATOM.{17}A\s{0,3}$each\b/ && $chain eq " ")
        {
        print PDBOUT "$line";
        }
     elsif($line =~ /^\s{0,4}ENDMDL/)
	{
	#print PDBOUT "END";
	last READLINE
	}
     elsif($line =~ /^\s{0,4}END/)
	{
	#print PDBOUT "END";
	last READLINE
	}
     else {}
  }
 }
print PDBOUT "END";
 close PDBOUT;

undef(@readpdb);
}
}

sub savehits #makes the results file containg all the hits
{

my @results=@{$_[0]};
my $writetype=$_[1];
my $resultsdir=$_[2];
my $runident=$_[3];
my $runcount=$_[4];
my $resultcnt=@results;

$runident=~s/:/_/g;
#print "0:$_[0],1:$_[1],2:$_[2],3:$_[3]\n";
my $filename=$resultsdir.$runident."\.smh";
print "$filename\n"; 
#my$filename="results:$seedpdb\_$iterations\_$dbfile";
#print "RES:$resultsdir\n";
#print "SAVEHITSFILENAME: $filename\n";
if ($writetype eq "a")
  {
  #print "savehits called append mode\n";
  open(OLDFILE,"<$filename");
  @oldresults=<OLDFILE>;

  #pop @oldresults;
  pop @oldresults;
  pop @oldresults;
  pop @oldresults;

  open(RESULTS,">$filename");
  for $line (@oldresults)
    {
    print RESULTS "$line";
    }
  close RESULTS;
  close OLDFILE;
 open(RESULTS,">>$filename");#append the file just created
  }

elsif ($writetype eq "b")
  {
  $filename=$resultsdir.$runident."-b"."\.smh";
  #print "savehits called";
  #print "file created: $filename\n";
open(RESULTS,">$filename"); # open new file
print RESULTS  "<XML>\n<smh type='bhits'>\n";
  }
elsif ($writetype eq "m")
  {
  $filename=$resultsdir.$runident."-m"."\.smh";
  #print "savehits called";
  #print "file created: $filename\n";
open(RESULTS,">$filename"); # open new file
print RESULTS  "<XML>\n<smh type='mhits'>\n";
  }
else
  {
  #print "savehits called";
  #print "file created: $filename\n";
open(RESULTS,">$filename"); # open new file
print RESULTS  "<XML>\n<smh type='hits'>\n";
  }




for (my $i=0;$i<$resultcnt; $i++)
  {
  my %each=%{$results[$i]};

  if(!$writetype eq "b")
  	{
	#print "ID:$each{file}, Binds:$each{metal} \n";
	}
  my@resnum=@{$each{resnum}};
  my@restype=@{$each{restype}};
  my@mresnum=@{$each{mresnum}};
  my@mrestype=@{$each{mrestype}};
 
  my $restotc=@{$each{resnum}};
  #my $hitid=$each{pdbid}.$each{resnum}[0]."\-".$each{resnum}[$restotc-1];
 
  
  print RESULTS "\n<hit id=$each{file} run=$each{pdbin}>\n";
  #print RESULTS "\n<hit id=$each{pdbid}>\n";
  print RESULTS "<iteration>".$each{iteration}."<\\iteration>\n";
  print RESULTS "<chain>".$each{chain}."<\\chain>\n";
if ($each{motif})
	{
  	print RESULTS "<motif>".$each{motif}."<\\motif>\n" ;
   	} 

if ($each{metal})
	{
	
	my @metal=@{$each{metal}};
	my $metcnt=@metal;
	#print "$metcnt METALs FOUND\n";
	for (my $j=0; $j<$metcnt;$j++)
		{
		my %metal=%{$metal[$j]};
		my $metal=$metal{met};
  		print RESULTS "<binds>".$metal."<\\binds>\n" ;
		}
   }
my $bvalue=$each{bvalue};
print RESULTS "<bindvalue>$bvalue<\\bindvalue>\n";  
#if ($each{scop}){
  #print RESULTS "<scop>".$each{scop}."<\\scop>\n" ;
  # }
for($r=0; $resnum[$r]; $r++)
    {
     print RESULTS "<res>"."$resnum[$r]"."<\\res>\n";
    }
if ($each{ligand})
   {	
	my @ligand=@{$each{ligand}};
	my $ligcnt=@ligand;
	for (my $i=0;$i<$ligcnt;$i++)
		{
		my %ligand=%{$ligand[$i]};
		print RESULTS "<ligand>\n"."<resid>".$ligand{resid}."<\\resid>\n"."<res>".$ligand{res}."<\\res>\n"."<\\ligand>\n" ;
   		}
   } 
    print RESULTS "<\\hit>\n";
  }
if ($writetype eq "a")
  {
 print RESULTS "<runcount>$runcount<\\runcount>\n";
	}
 print RESULTS  "\n<\\smh>\n<\\XML>\n";
 close RESULTS;
}
 
sub savefamilies
{
my $resultsdir=$_[1];
my $runident=$_[2];
 my$filename=$resultsdir.$runident."\.smf";
 
 #print "file created: $filename\n";
 open(RESULTS,">$filename"); # open new file

 print RESULTS  "<XML>\n<smh type='families'>\n";

 
my@results=@{$_[0]};
 for (my $i=0; $results[$i]; $i++)
  {
  my%each=%{$results[$i]};
  my@pdblist=@{$each{pdblist}};
  my$scop=$each{scop};
	my $bvalue=$each{bvalue};

 
  my $restotc=#my$filename="results:$seedpdb\_$iterations\_$dbfile";@{$each{resnum}};
 
  print RESULTS "\n<scopclass id=$scop>\n";
  print RESULTS "<bindvalue>$bvalue<\\bindvalue>\n";  
 my $pdbcount=@pdblist;
for($r=0; $r<$pdbcount; $r++)
    {
     print RESULTS "<id>"."$pdblist[$r]"."<\\id>\n";
    }

    print RESULTS "<\\scopclass>\n";
  }

 print RESULTS  "<\\XML>\n<\\smh>\n";
 close RESULTS;
}

sub addtofile
{
@filetoopen=$_[0];
@hitstoadd=@{$_[1]};

my$filename="/home/dwoodhead/wrkdir/results/"."hitsxml\-".$seedpdb."\-".$dbfile."\-".$rmscutoff."\-".$iterations."\.smh";
 #my$filename="results:$seedpdb\_$iterations\_$dbfile";

open(RESULTS,"<$filename");


}

sub saverms
{
my@hitlist=@_;

#$filename ="seqences\_$seedpdb\_$dbfile";
my$filename=$resultsdir.$runident."\.smr";
open(MATRIXFILE,">$filename"); # open new sequence file
my$hitlistcnt=@hitlist;

print MATRIXFILE  "<XML>\n<smh type='rms'>\n";
for (my$x=0; $x<$hitlistcnt; $x++)
{
	my%hit=%{$hitlist[$x]};
	
	my %rmslist=%{$hit{match}};
	@rmslist=%rmslist;
	#print "@rmslist\n";
	
	@rmslist=%rmslist;
	#print "rmslist array : \n";
	foreach (@rmslist) {
	
 	# print "$_\n";
			}
	
	
	$hitcnt=@rmslist;
	@checkarray;
	for ($i=0;$i<$hitcnt;$i++)
		{
		
		my@res=@{$hit{resnum}};
		$rescnt=@res;
		my$ID1=$hit{file};
		#print "ID1=$ID1\n";
		my$ID2=$rmslist[$i];
		#print "ID2=$ID2\n";
		$i++;
		my$rms=$rmslist[$i];
		#print "@rmslist\n";
		
		#$checkflag=0;
		#for $check(@checkarray)
			#{
			#if ($check eq "$ID1.$ID2")
			#	{$checkflag=1;}
			#else{}
			#}
		#if (!$checkflag)
		#{
		print MATRIXFILE "<RMS ID=".$ID1." ID=".$ID2 ."> ".$rms." <\\RMS>\n";
		#print "RMS saved: $ID1 , $ID2 , $rms\n";
		#my$checkid=$ID1.$ID2;
		#push (@checkarray,$checkid);
		#$checkid=$ID2.$ID1;
		#push (@checkarray,$checkid);
		#}
	}
}
print MATRIXFILE "<\\smh><\\XML>\n";
close MATRIXFILE;

}

sub saveproteins
{
my@pdblist=@{$_[0]};
my$dir=$_[1];

#$filename ="seqences\_$seedpdb\_$dbfile";
my$filename=$resultsdir."seq\-".$runident."\.pro";
open(SEQFILE,">$filename"); # open new sequence file

print SEQFILE "<xml>\n<smh type='protien'>\n";



$pdblist=@pdblist;

for ($i=0; $i<$pdblist; $i++)
  {
  my $seqdone=0;
  my $pdb=lc($pdblist[$i]{pdbid});
  my $seqnumber=@seqrecord;


  for (my $i=0; $i<$seqnumber;$i++)
    {
    $record=$seqrecord[$i];
    #print "$i record:$record:$pdb:pdb\n";
       if ($record eq $pdb)
          {
            #print "$pdb sequence already written\n";
            $seqdone=1;
          }
       else{}
       }
  if (not $seqdone)
   {
   #print "$pdb being written and added sequence to record\n";
            push(@seqrecord, $pdb);
            #print "@seqrecord\n";
   }
      }
for my$pdb (@seqrecord)
  {
  #print "pdb:$pdb\n";
  my$openpdb ="$dir$pdb.pdb";
  open(PDBFILE, $openpdb);
  @readpdb = <PDBFILE>;
  close PDBFILE;
  # Parse the record types of the PDB file
  my %recordtypes = parsePDBrecordtypes(@readpdb);

  # Extract the amino acid sequences of all chains in the protein
  my @chains = extractSEQRES( $recordtypes{'SEQRES'} );

  # Translate the 3-character codes to 1-character codes, and print
  
  
  print SEQFILE ">$pdb\n";
  foreach my $chain (@chains)
    {
    #print ">>>>$pdb\n";
    #print "****chain $chain **** \n\n";
    #print "$chain\n";
    print SEQFILE iub3to1($chain);
    }
   print SEQFILE "\n\n";
  }
 close SEQFILE;

 return;
 }

sub savedata
{

my $filename=$resultsdir."$runident.data";
my $filename2=$resultsdir."$runident AAType.txt";
my $filename3=$resultsdir."$runident ACC.txt";
my $filename4=$resultsdir."$runident ConRes.txt";
my $filename5=$resultsdir."$runident AAsize.txt";
my $filename6=$resultsdir."$runident SS.txt";


print ("file created: $filename\n");
open(RESULTS,">$filename"); # open new file
 
print ("file created: $filename2\n");
open(RESULTS2,">$filename2"); # open new file

 
print ("file created: $filename3\n");
open(RESULTS3,">$filename3"); # open new file

print ("file created: $filename4\n");
open(RESULTS4,">$filename4"); # open new file

print ("file created: $filename5\n");
open(RESULTS5,">$filename5"); # open new file

print ("file created: $filename6\n");
open(RESULTS6,">$filename6"); # open new file

print RESULTS  "<XML>\n<smh type='data'>\n";

 
my@results=@{$_[0]};
my $rescnt=@results;
for (my $i=0; $i<$rescnt; $i++)
	{
  	my%each=%{$results[$i]};
  	my @pdblist=@{$each{pdblist}};
  	my $scop=$each{scop};
	my $restotc=@{$each{resnum}};
 
  	print RESULTS "\n<scopclass id=$scop>\n";
  
 	my $pdbcount=@pdblist;
	for($r=0; $r<$pdbcount; $r++)
    		{
     		print RESULTS "<id>"."$pdblist[$r]"."<\\id>\n";
    		}
	my $alignfile=$each{alignfile};
    	print RESULTS "\t<Align>$alignfile<\\Align>\n";
    	print RESULTS "\t<AAtypedata>\n";
	print RESULTS2 "$scop\n";
	print RESULTS3 "$scop,";
	print RESULTS4 "$scop\n";
	print RESULTS5 "$scop\n\n";
	print RESULTS6 "$scop,";

	my @resdata=@{$each{resdata}};
	my $resdatacnt=@resdata;
		#print "resdata: @resdata\n";
		for (my $i=0;$i<$resdatacnt;$i++)
			{
			my %res=%{$resdata[$i]};
			my @res=%res;
			#print "RESDATA:@res\n";
			print RESULTS 	"\t\t<AA>\n";
			my $pos=$res{pos};
			my $res=$res{res};
			print RESULTS "\t\t\t<Pos>$pos<\\Pos>\n";
			print RESULTS "\t\t\t<Res>$res<\\Res>\n";
			my $L=$res{'perL'};
			my $S=$res{'perS'};
                    	my $Ne=$res{'perNe'};
			my $Pl=$res{'perPl'};
			my $P=$res{'perP'};
			my $N=$res{'perN'};
			my $G=$res{'perG'};
			print RESULTS2 "$pos,$Ne,$Pl,$P,$N,$G,\n";
			print RESULTS "\t\t\t<size>\n";
			if($L){print RESULTS "\t\t\t\t<L>$L<\\L>\n";}
			if($S){print RESULTS "\t\t\t\t<S>$S<\\S>\n";}
			print RESULTS5 "$L,$S,\n";
			print RESULTS "\t\t\t<\\size>\n";
			print RESULTS "\t\t\t<group>\n";
			if($Ne){print RESULTS "\t\t\t\t<Nve>$Ne<\\Nve>\n";}
			if($Pl){print RESULTS "\t\t\t\t<Pve>$Pl<\\Pve>\n";}
			if($P){print RESULTS "\t\t\t\t<P>$P<\\P>\n";}
			if($N){print RESULTS "\t\t\t\t<N>$N<\\N>\n";}
			if($G){print RESULTS "\t\t\t\t<G>$G<\\G>\n";}
			print RESULTS "\t\t\t<\\group>\n";
				
			print RESULTS "\t\t<\\AA>\n";
			}
	
	print RESULTS "\t<AAtypedata>\n";
 	print RESULTS "\t<SSI>\n";
	print RESULTS "\t\t<SSfile>$each{SSfile}<\\SSfile>\n";
	my %ssdata=%{$each{SSdata}};
	my $upl=$ssdata{uplength};
	my $upd=$ssdata{updist};
	my $dwnl=$ssdata{downlength};
	my $dwnd=$ssdata{downdist};
	my $sss=$ssdata{ssscore};
	my $ssline1=$ssdata{ssseq1};
	my $ssline2=$ssdata{ssseq2};
	print RESULTS "\t\t<USSS>\n\t\t\t<Length>$upl<\\Length>\n\t\t\t<Distance>$upd<\\Distance>\n\t\t<\\USSS>\n";
	print RESULTS "\t\t<MSS>\n\t\t\t<Score>$sss<\\Score>\n\t\t<\\MSS>\n";
	print RESULTS "\t\t<DSSS>\n\t\t\t<Length>$dwnl<\\Length>\n\t\t\t<Distance>$dwnd<\\Distance>\n\t\t<\\DSSS>\n";
    	print RESULTS "\t<\\SSI>\n";
	print RESULTS6 "$upl,$upd,$sss,$dwnl,$dwnd\n";
	print RESULTS6 "$ssline1\n";
	print RESULTS6 "$ssline2\n";
	#my $solacc=$each{solacc};
	#my $acc=$each{acc};
	print RESULTS "\t<SA>\n";
	#print RESULTS "\t\t<SACCS>$acc<\\SAACS>\n";
	#my@acc=split(//,$acc);
	#for my $AR(@acc)
	#{
	#$AR=altaccscore($AR);
	#print RESULTS3 "$AR,";
   	#}
	#print RESULTS3 "\n";
	#print RESULTS "\t\t<SACC>$solacc<\\SACC>\n";
	my $res3p;
	my $res3a=" ,";
	my $res3c=" ,";
	@sable=@{$each{sable}};
	$sablecnt=@sable;
	for (my $i=0; $i<$sablecnt; $i++)
		{
		%sableres=%{$sable[$i]};
		print RESULTS 	"\t\t<SABRES>\n";
			my $pos=$sableres{pos};
			my $res=$sableres{res};
			my $sacc=$sableres{acc};
			my $conf=$sableres{conf};
			print RESULTS "\t\t\t<Pos>$pos<\\Pos>\n";
			$res3p=$res3p.$pos.",";
			print RESULTS "\t\t\t<RES>$res<\\RES>\n";
			print RESULTS "\t\t\t<ACC>$sacc<\\ACC>\n";
			$res3a=$res3a.$sacc.",";
			print RESULTS "\t\t\t<CONF>$conf<\\CONF>\n";
			$res3c=$res3c.$conf.",";
		print RESULTS 	"\t\t<\\SABRES>\n";
		}
	print RESULTS "\t<\\SA>\n";
	print RESULTS3 "$res3p\n$res3a\n$res3c\n";
	print RESULTS "\t<CONDE>\n";
	my @conserved=@{$each{conserved}};
	my $concnt=@conserved;
	my $Dflag=0;
	my $Eflag=0;
	for (my $i=0; $concnt>$i; $i++)
	{
	my %con=%{$conserved[$i]};
	#my @con=%con;
	#print "ConData:@con\n";
	my $dist=$con{dist};
	my $perD=$con{percentD};
	my $perE=$con{percentE};
	print RESULTS "\t\t<CON>\n\t\t\t<dist>$dist<\\dist>\n\t\t\t<D>$perD<\\D>\n\t\t\t<E>$perE<\\E>\n\t\t<\\CON>\n";
	
	
	}
	my $groupnumber=11;
	my $groupsize=10;
	my $groupbttm=1;
	my $grouptop=$groupsize;
	my $groupline;
	my $groupcntline;
	
	for (my $i=0;$i<$groupnumber;$i++)
		{
		my $cnt;
		if ($i==($groupnumber-1))
			{
			$cnt=testconserved(\@conserved,'1000',$groupbttm);
			$groupline=$groupline."$groupbttm +,";
			}
		else 
			{
		$cnt=testconserved(\@conserved,$grouptop,$groupbttm);
		$groupline=$groupline."$groupbttm - $grouptop,";
			}
		$groupcntline=$groupcntline."$cnt,";
		$groupbttm=$groupbttm+$groupsize;
		$grouptop=$grouptop+$groupsize;
		}
	print RESULTS4 "$groupline\n$groupcntline\n";
	print RESULTS "\t<\\CONDE>\n";
    	print RESULTS "<\\scopclass>\n";
  	}
 print RESULTS  "<\\XML>\n<\\smh>\n";
 close RESULTS;
}

sub addtofile
{
@filetoopen=$_[0];
@hitstoadd=@{$_[1]};

my$filename="/home/dwoodhead/wrkdir/results/"."hitsxml\-".$seedpdb."\-".$dbfile."\-".$rmscutoff."\-".$iterations."\.smh";
 #my$filename="results:$seedpdb\_$iterations\_$dbfile";

open(RESULTS,"<$filename");


}

1;