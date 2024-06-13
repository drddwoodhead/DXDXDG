#use SOAP::Lite;
use IO::File; 


sub addscopclass
{
my @hitlist=@{$_[0]};
my $flag=$_[1];
my $scopdb=$_[2];
my $hitlist=@hitlist;
my @outlist;
#print "Classifying $hitlist hits\n|0%......................|50%.....................|100%\n";
my $percentage=int($hitlist/50);
#print "percentage:$percentage\n";
open(SCOPFILE, $scopdb);
my @readdb = <SCOPFILE>;
my %dbhash=scopsearchhash(\@readdb);
while ( ($k,$v) = each %dbhash ) {
   # print "$k => $v\n";
}
close SCOPFILE;
undef @readdb;
for  (my $i=0 ; $i<$hitlist; $i++)
	{
#	if ($i % $percentage==0)
#	{
#	print "*";
#	}
	my%each=%{$hitlist[$i]};
	my $pdbid=$each{pdbid};
	my $chain=$each{chain};
        my @res=@{$each{resnum}};
        my $res1=$res[0];
        my $cnt=@res;
        my $res2=$res[$cnt-1];
  #	print "finding scop for :$pdbid:$chain:$res1:\-:$res2:\($i of $hitlist\)\n";
	my $class=findclass($pdbid,$chain,$res1,$res2,$scopdb,\%dbhash);
#	print "result:$class\n";
	if ($class)
	  {
	  my $source=$scopdb;
	  $source=~s/^.*\///;
	  $each{predicted}=0;
	  $each{scopsource}=$source;
	  $each{scop}=$class;
	  }
	elsif ($class eq "not found" && $flag ne 'nossm')
	  {
	  %ssmresult=predictclass($pdbid,$chain,$res,$res2,$scopdb);
	  $resultid=$ssmresult{targetid}.":".$result{targetchain};
	  $resultsun=$ssmresult{sunid};
	  $resultscop=$ssmresult{scop};
          $qscore=$ssmresult{Q};
          @ssepairs=@{$resultssm{ssepairs}};
          %nearestresult=%{$ssepairs[0]};
          $ssedistance=$nearestresult{distance};
          $ssetype=$nearestresult{type};
	  my $source="predicted by similarity to $resultid by ssm with Q score of $qscore";
	  $each{predicted}=1;
	  $each{scopsource}=$source;
	  $each{ssmmatch}=\%ssmresult;
	  $each{scop}=$class;
	  }
	push(@outlist,\%each);
	}
	print "\n";
	return @outlist;
}

sub findclass #get the scop class from the pdb , chain and residues or from the sun id.

{
my @args=@_;
my $argcnt=@args;
my $id;
if ($argcnt>2)
{
#print "$argcnt args found\n";
my $pdbid=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];
my $scopdb=$_[4];
my %dbhash=%{$_[5]};
$id=findclassid($pdbid,$chain,$res1,$res2,$scopdb,\%dbhash);
}
else
{
$id=$_[0];
$scopdb=$_[1];
}
#print "using sunid : $id\n";
open(SCOPFILE, $scopdb);
my @readdb = <SCOPFILE>;

for $line(@readdb)
  {  
if( $line =~ /.*(\w\.\d{1,3}\.\d{1,3}\.\d{1,3})\s$id.*/)
  	{
	close SCOPFILE;
	my $scop=$1;
        return $scop;
        }
  }
close SCOPFILE;
return "not found";
}

sub findclassid #looks for pdbid in the scop database and returns its classification.
 #arguments are the 4 character pdb id
{
my $pdbid=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];
my $scopdb=$_[4];
my %dbhash=%{$_[5]};
 ## arguments

#print "arguments :$pdbid, :$chain, :$res1-$res2\n";
#my $scopdb="/home/dwoodhead/scop/dir.cla.scop.txt_masterfile-pre-scop-2007-06-01.out";
if(!open(SCOPFILE ,$scopdb))
{print "SCOP database file not found\n$scopdb\n";
}
else {#print "SCOP database loaded\n $scopdb\n";
}
open(SCOPFILE, $scopdb);
my @readdb = <SCOPFILE>;


my $linecnt=@readdb;

my $index=$dbhash{$pdbid};
#print "$pdbid found index:$index\n";
if($index!="")
{
for (my $i=$index,$i<$linecnt,$i++)
  { 
	$line=$readdb[$i];
  if( $line =~ /.{8}$pdbid\s$chain:\d/)
  	{
	#print "ID Chain:range\n";
	#print "pdbfound in scop database (chain/res)\n";
	#print "$line";
  	$line=~m/.{8}$pdbid\s$chain:(\d{1,4})-(\d{1,4})/;
  	my$low=$1;
  	my$high=$2;
	#print "low:$low res:$res1-$res2 high:$high\n";
  	if($res1>=$low&&$res2<=$high)
    		{$line=~m/.*\w\.\d{1,3}\.\d{1,3}\.\d{1,3}\s(\d{1,6}).*/;
	   	my$scopclass=$1;
		#print "FOUND:$scopclass\n";
		undef @readdb;
		close SCOPFILE;
    		 return $scopclass;
   		}  
	else{}
    	}
  elsif( $line =~ /.{8}$pdbid\s{1}\d{1,4}-\d{1,4}/)
    	{
	#print "ID range\n";
	#print "$line";
	#print "pdbfound in scop database (res)\n";
	$line=~m/.{8}$pdbid\s{1}(\d{1,4})-(\d{1,4})/;
 	my $low=$1;
  	my $high=$2;
	#print "low:$low res:$res1-$res2 high:$high\n";
     	 if($res1>=$low&&$res2<=$high)
    		{$line=~m/.*\w\.\d{1,3}\.\d{1,3}\.\d{1,3}\s(\d{1,6}).*/;
	 	 my$scopclass=$1;
		undef @readdb;
		close SCOPFILE;
	   	return $scopclass;
    		}
   	else{}
    	}
  elsif( $line =~ /.{8}$pdbid\s{1}$chain/)
  	{
	#print "ID Chain\n";
	#print "pdbfound in scop database (chain)\n";
	#print "$line";
	#my$scopclass=substr($line,16,8);
	$line=~m/.*\w\.\d{1,3}\.\d{1,3}\.\d{1,3}\s(\d{1,6}).*/;
	my$scopclass=$1;
	#print"$line\n";
	#print "scop class found:$scopclass\n";
	undef @readdb;
	close SCOPFILE;
	return $scopclass; 
	}
  else
	{
	#print ".";
	}
  }

if ($chain=="-")
{
for $line(@readdb)
  { 
  if($line =~ /.{8}$pdbid\s\w:{1}\d{1,4}-\d{1,4}/)
    	{
	#print "ID range\n";
	#print "$line";
	#print "pdbfound in scop database (res)\n";
	$line=~m/.{8}$pdbid\s{1}(\d{1,4})-(\d{1,4})/;
 	my $low=$1;
  	my $high=$2;
	#print "low:$low res:$res1-$res2 high:$high\n";
     	 if($res1>=$low&&$res2<=$high)
    		{$line=~m/.*\w\.\d{1,3}\.\d{1,3}\.\d{1,3}\s(\d{1,6}).*/;
	 	 my$scopclass=$1;
		undef @readdb;
		close SCOPFILE;
	   	return $scopclass;
    		}
   	else{}
    	}

}

for $line(@readdb)
  { 
  if($line =~m/.{8}$pdbid\s{1}\w:\d/)
  	{
	#print "ID Chain\n";
	#print "pdbfound in scop database (chain)\n";
	#print "$line";
	#my$scopclass=substr($line,16,8);
	$line=~m/.*\w\.\d{1,3}\.\d{1,3}\.\d{1,3}\s(\d{1,6}).*/;
	my$scopclass=$1;
	#print"$line\n";
	#print "scop class found:$scopclass\n";
	undef @readdb;
	close SCOPFILE;
	return $scopclass; 
	}
  else{}
  }
}
}
undef @readdb;
close SCOPFILE;
return "not found"; 
}


sub predictclass
{
my $pdbid=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];


my @xmlresult=runssm($pdbid,$chain,$res1,$res2);
my @ssmhits=parsessmxml(\@xmlresult);
my %match=pickresult(\@ssmhits);
return %match;

}

sub runssm
{
my $service = SOAP::Lite
  ->uri('urn:msd_soap_service')
  ->proxy('http://www.ebi.ac.uk/msd-srv/msdsoap'); # Or:
# -> service('http://www.ebi.ac.uk/msd-srv/docs/api/msd_soap_service.wsdl');

my $pdbinput=$_[0];
my $chain=$_[1];
my $firstres=$_[2];
my $lastres=$_[3];

my $rangeselect=50;
my $range1=$firstres-$rangeselect;
if ($range1<0)
{$range1=0;}
else{}
my $range2=$lastres+$rangeselect;

my $chainselect=$chain.":".$chain.$range1."-".$range2;


 my $query = '<SSMInput> <query> <type>PDB entry</type> <pdbcode>'.$pdbinput.'</pdbcode>'.
'</query> <target> <type>PDB archive</type>'.
'</target> <selection> <type>Chain(s)</type>'.
'<chains>'.$chainselect.'</chains> </selection> <percent1>70</percent1>'.
'<percent2>70</percent2> <sepchains>Yes</sepchains>'.
'<connectivity>Yes</connectivity> <bestmatch>Yes</bestmatch>'.
'<uniquematch>Yes</uniquematch> <precision>Normal</precision>'.
'<sorting>RMSD</sorting> </SSMInput>';

my $sessionid = 'mysessionnumber';
my $array = [$query , $sessionid];
print "\nCalling method: msdSSM()\n";
my $result1 = $service->msdSSM(SOAP::Data->name('numofpars' => 2),
   SOAP::Data->name('inparams' => $array));

unless ($result1->fault) {
  print "\nmsdSSM() called successfully.\n";
  #print $result1->result();
}else{
  print join ', ',
  $result1->faultcode,
  $result1->faultstring;
}
for (my $i=1; $i > 0 ; $i++)
  { my $result2 = $service->msdGetSoapBase64Data(SOAP::Data->name('diyf')->type('xsd:int')->value(1),
   SOAP::Data->name('sessionid')->type('xsd:string')->value($sessionid));
   unless ($result2->fault || ($result2->faultstring != 'Result file not ready')) {
   $i=-1;
   #print $result2->result();
   @xmlarray=$result2->result();
   $fh = new IO::File "> $resultsdir.$pdbinput\:$chain\:$range1\-$range2.xml";
   if (defined $fh){
    print $fh $result2->result();
    $fh->close;
    print "\nResult file created successfully:\n$fh\n";
   }
   }else{
   print join ', ',
   "\nXML File not ready...retry: $i ",
   $result2->faultcode,
   $result2->faultstring;
   }
   sleep(5);
}
   my $result3 = $service->msdSSMPurge(SOAP::Data->name('sessionid')->type('xsd:string')->value($sessionid));
   unless ($result3->fault) {
   print "\nmsdSSMPurge() called successfully.\n";
   print $result3->result();
  }else{
   print join ', ',
   $result3->faultcode,
   $result3->faultstring;
 }
#print @xmlarray;
return @xmlarray;
}

sub parsexmlssm   #parses the ssm xml file
{
my$openfile=$_[0];
print "locatoin of file :$openfile\n";
open(FILE, $openfile);
@readfile = <FILE>;
close FILE;
print $readfile[0];
 if ($readfile[0] eq "<SSMResults>\n")
 {
    print "xml type = hits\n";
   @outarray=parsessmxml(\@readfile);
return @outarray;
 }
   
 else
 {
   print "no xml file found please check format and try again\n";
 }
}

sub parsessmxml
{
my @readfile=@{$_[0]};

my @outarray;

my $inmatch='0';



#print "parse ssm xml initialised\n";

$linecount=@readfile;
#print "line count: $linecount\n";
  my @targetsunids;
  my @querysunids;
  my $Q;
  my $P;
  my $Z;
  my $targetid;
  my $targetchain; 
  my @ssepairlist;
  my %sse;
  my %tsse;
  my %qsse;

for (my $i=0;$linecount>$i;$i++)
{
$line=$readfile[$i];
 
 
  #print "$i not in hit\n";
if ($line =~m/<Match>\n/)
    #if ($line =~m/<hit id\=\w{4}\d{1,4}-\d{1,4}>\n/)
    {
    $Q=0;
    $P=0;
    $Z=0;
    $inmatch=1;
    #print "hit found:\n";
    #print "$line";
    }
	

 
MATCHLOOP: while($inmatch)
  {
  $line=$readfile[$i];
        		
  if ($line =~m/^\s*<Q-score>(.*)<\/Q-score>\n/)
	{
	$Q=$1;
	#print "q:$Q\n";
        }
  elsif ($line =~m/^\s*<P-score>(.*)<\/P-score>\n/)
        {
        $P=$1;
	#print "p:$P\n";
        }
  elsif ($line =~m/^\s*<Z-score>(.*)<\/Z-score>\n/)
        {
        $Z=$1;
	#print "z:$Z\n";
        }
  elsif ($line =~m/<Query>/)
        {
  	undef(@querysunids);
	#print "$i$line";
	$inquery=1;
        }
  elsif ($line =~m/<Target>/)
        {
	$targetid="";
  	$targetchain="";   
	undef(@targetsunids);
 	#print "$i$line";
        $intarget=1;
	}
  elsif ($line =~m/<RTMatrix>/)
        {
	#print "$i$line";
        $inrtmatrix=1;
        }
  elsif ($line =~m/<SSEMatch>/)
        {
	#print "$i$line";
        $inssematch=1;
          }
  elsif ($line =~m/^\s*<\/Match>/)
        {
        #print "end of hit found\n\n";
        my $pdbid=substr($hitid,0,4);
	#print "$targetid, $targetchain, Q=$Q,\n";
	#print "target sun ids:\n";
	for my $every(@targetsunids)
		  {
		  #print ":$every\n";
		  }
	my @tsunidref=@targetsunids;
	#print "query sun ids:\n";
	for my $every(@querysunids)
		{		
		#print ":$every\n";
		}
	my $ssepaircount=@ssepairlist;
	#print "$ssepaircount pairs found\n";
	my @qsunidref=@querysunids;
	my @ssepairref=@ssepairlist;
	undef(@ssepairlist);
        my%hit=("targetid",$targetid,"targetchain",$targetchain,"Q",$Q,"Z",$Z,"P",$P, "ssepairs",\@ssepairref,"targetsunids",\@tsunidref,"querysunids",\@qsunidref);
        push(@outarray,\%hit);
        $inmatch=0;
        }
        
	while($inquery)
		{
		#print "while loop active\n";
		$line=$readfile[$i];
        	#print "$i is in hit\n";
		if ($line =~m/<\/Query>/)
        		{
			#print "$line";
        		$inquery=0;
			last MATCHLOOP;
          		}
                elsif ($line =~m/^\s*<scopdomains>(.*)<\/scopdomains>/)
                			{
			#$numberofscopdomains=$1;
			}
		elsif ($line =~m/^\s*<sunid>(.*)<\/sunid>/)
			{
			$sunid=$1;
			#print "target sunid:$sunid\n";
			push(@querysunids,$sunid);
			}
		$i++;
		
		}
	while($intarget)
		{
		#print "while loop active\n";
		$line=$readfile[$i];
        	#print "$i is in hit\n";
		if ($line =~m/<\/Target>/)
        	        {
			#print "$i$line";
        		$intarget=0;
			last MATCHLOOP;
          		}
		elsif ($line =~m/^\s*<scopdomains>(.*)<\/scopdomains>/)
                			{
			#$numberofscopdomains=$1;
			}
		elsif ($line =~m/^\s*<sunid>(.*)<\/sunid>/)
			{
			$sunid=$1;
			#print "query sunid:$sunid\n";
			push(@targetsunids,$sunid);
			}
		elsif ($line =~m/^\s*<name>(.*)<\/name>/)
			{
			$targetid=$1;
			}
		elsif ($line =~m/^\s*<chainID>(.*)<\/chainID>/)
			{
			$targetchain=$1;
			}
      
		$i++;
		}
	while($inrtmatrix)
		{
		$line=$readfile[$i];
		if ($line =~m/<\/RTMatrix>/)
        		{
        		$inrtmatrix=0;
			#print "$line";
			last MATCHLOOP;
          		}
		$i++;
		}
SSEMLOOP: while($inssematch)
		{
		$line=$readfile[$i];
		if ($line =~m/<\/SSEMatch>/)
        		{
        		$inssematch=0;
			#print "$i$line";
			last MATCHLOOP;
          		}
		elsif ($line =~m/<SSEpair>/)
			{
			#print "$i$line";
            		$inssepair=1;
			}

			SSEPLOOP: while ($inssepair)
				{
				$line=$readfile[$i];
				if ($line =~m/<\/SSEpair>/)
					{
					$inssepair=0;
					#print "$i$line";
					my $ssecnt=@ssepairlist;
					#print"added another pair now $ssecnt\n";
					my %ssepush=%sse;
					push (@ssepairlist,\%ssepush);
					last SSEMLOOP;
					} 
				elsif($line =~m/<type>(.*)<\/type>/)
					{
					#print "$i$line";;
					$sse{type}=$1;
					}
				elsif($line =~m/<Query>/)
					{
					#print "Q$i$line";
            				$inquery=1;
					}
				elsif($line =~m/<Target>/)
					{
					#print "T$i$line";
            				$intarget=1;
					}
				SSEQLOOP: while ($inquery)
					{
					$line=$readfile[$i];
					if ($line =~m/<\/Query>/)
						{
						$inquery=0;
						#print "Q$i$line";
						my %qref=%qsse;
						$sse{query}=\%qref;
						last SSEPLOOP;
						}
					elsif($line =~m/<chainID>(.*)<\/chainID>/)
						{
						#print "q$i$line";
						$qsse{chain}=$1;
						}
          				elsif($line =~m/<length>(.*)<\/length>/)
						{
						#print "q$i$line";
						$qsse{length}=$1;
						}
					elsif($line =~m/<initRes>/)
						{
						#print "q$i$line";
						$intres=1;
						}
					elsif($line =~m/<endRes>/)
						{
						#print "q$i$line";
						$endres=1;
						}
					else
						{
						#print "_$i$line";
						}
					while ($intres)
						{
						$line=$readfile[$i];
						if ($line =~m/<\/initRes>/)
							{
							$intres=0;
							#print "E$i$line";
							last SSEQLOOP;
							}
						elsif($line =~m/<seqNum>(.*)<\/seqNum>/)
							{
							#print "q$i$line";
							my $res=$1;
							$qsse{res}=$res;
							#print"query residue:$qsse{res}\n";
							}
						elsif($line =~m/.*/)
							{
							#print "A$i$line";
							}
						$i++;
						}
					while ($endres)
						{
						$line=$readfile[$i];
						if ($line =~m/<\/endRes>/)
							{
							$endres=0;
							#print "E$i$line";
							last SSEQLOOP;
							}
						else
							{
							#print "L$i$line";
							}
						$i++;
						}

					$i++;
					}
								
				SSETLOOP: while ($intarget)
					{
					$line=$readfile[$i];
					if ($line =~m/<\/Target>/)
						{
						$intarget=0;
						#print "T$i$line";
						my %tref=%tsse;
						$sse{target}=\%tref;
						last SSEPLOOP;
						}
					elsif($line =~m/<chainID>(.*)<\/chainID>/)
						{
						#print "t$i$line";
						$tsse{chain}=$1;
						}
          				elsif($line =~m/<length>(.*)<\/length>/)
						{
						#print "t$i$line";
						$tsse{length}=$1;
						}
					elsif($line =~m/<initRes>/)
						{
						#print "t$i$line";
						$intres=1;
						}
					elsif($line =~m/<endRes>/)
						{
						#print "t$i$line";
						$endres=1;
						}
					while ($intres)
						{
						$line=$readfile[$i];
						if ($line =~m/<\/initRes>/)
							{
							$intres=0;
							#print "t$i$line";
							last SSETLOOP;
							}
						elsif($line =~m/<seqNum>(.*)<\/seqNum>/)
							{
							my $res=$1;
							$tsse{res}=$res;
							#print"targetresidue:$tsse{res}\n";
							}
						$i++;
						}
					while ($endres)
						{
						$line=$readfile[$i];
						if ($line =~m/<\/endRes>/)
							{
							$endres=0;
							#print "t$i$line";
							last SSETLOOP;
							}
						$i++;
						}
					$i++;
					}

				
				$i++;
				}
			$i++;
			}	
        	$i++;
        }
  
}
$outcount=@outarray;

#print "$outcount Matches found\n";

return @outarray;
}

sub sortssebynearest
{
my $pdbid=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];
my @matchlist=@{$_[4]};

my %sse;
my $distance;
my %match;
my @outlist;
my @sselist;
my @out;

my $matchcount=@matchlist;

for (my$i=0;$i<$matchcount;$i++)
{ 
%match=%{$matchlist[$i]};
@sselist=@{$match{ssepairs}};
undef(@out);

#print "Sorting SSEs by nearest to motif...\n Args :$pdbid, :$chain, :$res1-$res2\n";

my $ssecount=@sselist;

for (my $i=0; $i<$ssecount; $i++)
	{
	%sse=%{$sselist[$i]};
	my %query=%{$sse{query}};
	my $ssechain=$query{chain};
	my $ssestart=$query{res};
	my $sseend=$query{res}+$query{length};
	
	my %target=%{$sse{target}};
	my $tssechain=$target{chain};
	my $tssestart=$target{res};
	my $tsseend=$target{res}+$query{length};
	
	#print "Qsse $i of $ssecount, $ssechain:$ssestart-$sseend\n";
	#print "Tsse $i of $ssecount, $tssechain:$tssestart-$tsseend\n";
	if ($chain eq $ssechain)
		{
		 my $sub1=$res1-$ssestart;
		 my $sub2=$res1-$sseend;
		 my $sub3=$res2-$ssestart;
		 my $sub4=$res2-$sseend;
		 $sub1=sqrt($sub1*$sub1); 
		 $sub2=sqrt($sub2*$sub2); 
		 $sub3=sqrt($sub3*$sub3); 
		 $sub4=sqrt($sub4*$sub4); 
		 #print "$sub1, $sub2, $sub3, $sub4\n";
		 $distance=smallest($sub1,$sub2,$sub3,$sub4);
		#print "calculated distance is:$sub4istance\n";
		
		}
	else {
		$distance='999999';
		#print"different chains distance set to 999999\n"; 
		}
	my %sseref=%sse;
	#print "adding distance :$distance to sse\n";
	my $distref=$distance;
	$sseref{distance}=$distref;
	push (@out, \%sseref);
	}
@out= sort {$$a{distance} <=> $$b{distance}} @out;
my @sseout=@out;
$match{ssepairs}=\@sseout;
my %outmatch=%match;
push (@outlist,\%outmatch);
}
return @outlist;
}

sub ssmmatchfindsunid
{
my $pdbid=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];
my %match=%{$_[4]};

my @sselist=@{$match{ssepairs}};

my $sselistcount=@sselist;
#@sselist=sortssebynearest($pdbid,$chain,$res1,$res2,\@sselist);

#print "$sselistcount structural elements \n"; 
for (my $i=0; $i<$sselistcount; $i++)
{
my %nearestsse=%{$sselist[$i]};
my %target=%{$nearestsse{target}};

my $ssestart=$target{res};
my $sseend=$target{res}+$target{length};
my $matchid=$match{targetid};
my $chain=$match{targetchain};
my $type=$target{type};
#print "looking for sunid of$matchid,$chain,$ssestart-$sseend\n";
my $id=findclassid($matchid,$chain,$ssestart,$sseend);
#print "id found: $id\n";

if ($id && $id ne "not found")
{
#print"sunid:$id\n";
return $id;
}
else 
{
}
}
#print "sunid: not found\n";
return "not found";


}

sub pickresult
{
my @resultlist=@{$_[0]};
my $pdbid=$_[1];
my $chain=$_[2];
my $res1=$_[3];
my $res2=$_[4];

my $listcnt=@resultlist;

my @picklist;


for (my $i; $i<$listcnt; $i++)
{
my %match=%{$resultlist[$i]};

my $sunid=ssmmatchfindsunid($pdbid,$chain,$res1,$res2,\%match);
my @sselist=@{$match{ssepairs}};
my %nearestsse=%{$sselist[0]};
my $distance=$nearestsse{distance};
#print"sunid from ssmmatch, :$sunid, (sse used $distance residues)\n";
$match{sunid}=$sunid;
my $scopclass=findclass($match{sunid});
$match{scop}=$scopclass;
#print "scop class found:$scopclass\n";
if ($match{scop} && $match{scop} ne "not found")
{
#print "added to list\n";
push (@picklist, \%match);
}
else {}
}
$pickcount=@picklist;
#print "$pickcount possible matches\n";
if ($pickcount>0)
{
return %{$picklist[0]};
}

else 
{
%notfound=("scop","not found");
return %notfound;
}
}

sub noscoplist
{
@familydata=@{$_[0]};

$familycount=@familydata;
my @pdbids;

for (my $i=0; $i<$familycount; $i++)
{
my %family=%{$familydata[$i]};
if ($family{scop} eq "not found")
{
@pdbids=@{$family{pdblist}};
}
else
{}

}
return @pdbids;
}

sub sortbyclass
{
my @classlist;

my @hitlist=@{$_[0]};
$hitlist=@hitlist;

for  (my $i=0 ; $i<$hitlist; $i++)
	{
	my%each=%{$hitlist[$i]};
	my $pdbid=$each{file};
	my $bvalue=$each{bvalue};
	if ($each{predicted}==1)
		{
		$pdbid=$pdbid."(p)";
		}
	my $scop=$each{scop};
	if ($scop eq "")
	{$scop="not found";}
	my $classflag=0;
	my %selectedclass=0;
	my $arrayref=0;

my$classlist=@classlist;	

	for  ($j=0 ; $j<$classlist; $j++)
	{
	my%classentry=%{$classlist[$j]};
	my$recordedscop=$classentry{scop};
	#print"i=$i j=$j comparison $scop == $recordedscop\n";
	if($scop eq $recordedscop)
	 {
          #print "TRUE setting flag to 1\n";
	  $classflag=1;
	  %selectedclass=%classentry;
	  $arrayref=$j;	
	 }
	else {
		#print "flag is $classflag\n";
		};
	}
if ($classflag)
{
my@pdblist=@{$selectedclass{pdblist}};
push(@pdblist,$pdbid);
$selectedclass{pdblist}=\@pdblist;
$classlist[$arrayref]=\%selectedclass;
#print "class true adding $pdbid to class $selectedclass{scop}\n";
}
else
{
my%newclass;
$newclass{scop}=$scop;
$newclass{bvalue}=$bvalue;
my@pdblist;
push(@pdblist,$pdbid);
$newclass{pdblist}=\@pdblist;
#$newclass{class}=;#placeholder for the written classification.
#print "new class found, creating class $scop, with pdb $pdbid\n"; 

push(@classlist,\%newclass)
}

}	

@classlistsrt = sort {$$a{scop} cmp $$b{scop}} @classlist;

return @classlistsrt;
}


sub saveformated
{
my@scoplist=@{$_[0]};
my @bindinglist=@{$_[1]};
my @hitlist=@{$_[2]};

my $hitlistcnt=@hitlist;

 my$filename=$resultsdir."formattedxml\-".$runident."\.smh";
 
# print ("file created: $filename\n");
 open(RESULTS,">$filename"); # open new file

 print RESULTS  "<XML>\n<smh type='families'>\n";

 

 for (my $i=0; $scoplist[$i]; $i++)
  {
  my%scopclass=%{$scoplist[$i]};
  my@scoppdblist=@{$scopclass{pdblist}};
  my$scop=$scopclass{scop};

 
  my $restotc=#my$filename="results:$seedpdb\_$iterations\_$dbfile";@{$each{resnum}};
 
  print RESULTS "\n<scopclass id=$scop>\n";
  

 for (my $j=0; $bindinglist[$j]; $j++)
  {
  my%metalbound=%{$bindinglist[$j]};
  my@bindingpdblist=@{$metalbound{pdblist}};
  my$metal=$metalbound{metal};

  my@idlist;

  my $scoppdbcount=@scoppdblist;
  for(my $r=0; $r<$scoppdbcount; $r++)
    	{
	my $pdbidscop=@scoppdblist[$r];
     
        my $bindingpdbcount=@bindingpdblist;
	for(my $s=0; $s<$bindingpdbcount; $s++)
		{
		my $pdbidbound=@bindingpdblist[$s];
		if ($pdbidscop eq $pdbidbound)
			{
			push ( @idlist, $pdbidbound) ;
			}
		else {}
		}
    	}
    my $idcount=@idlist;

    if ($idcount>0)
   {
    print RESULTS "\n<metalbound id=$metal>\n";
    for $each(@idlist)
	{
print RESULTS "<hit>";
	print RESULTS "<id>"."$each"."<\\id>";
	 for (my $i=0 ; $i<$hitlistcnt ; $i++)
		{
		my %every=%{$hitlist[$i]};
		my $everyid=$every{file};
		#print RESULTS "$everyid==$each";
		if ($everyid eq $each)
			{
			my $motif=$every{motif};
			print RESULTS "<motif>"."$motif"."<\\motif>";
			}
		}

	print RESULTS "<\\hit>\n";
	}
    print RESULTS "<\\metalbound>\n";
   }
  }
 print RESULTS "\n<\\scopclass>\n";
}
 print RESULTS  "<\\XML>\n<\\smh>\n";
 close RESULTS;
}

sub sortbybinding1
{
#print "sort by binding\n";
my @bindinglist;

my @hitlist=@{$_[0]};
$hitlist=@hitlist;

for  (my $x=0 ; $x<$hitlist; $x++)
	{
	my %each=%{$hitlist[$x]};
	my $pdbid=$each{file};
	my @metal=@{$each{metal}};
	my $bindingflag=0;
	my %selectedbinding=0;
	my $arrayref=0;

my$bindinglist=@bindinglist;	

	for  (my $j=0 ; $j<$bindinglist; $j++)
	{
	my%bindingentry=%{$bindinglist[$j]};
	my$recordedmetal=$bindingentry{metal};
	#print"i=$i j=$j comparison $metal == $recordedmetal\n";
	
	for $metal(@metal)
	{
	if($metal eq $recordedmetal)
	 {
          #print "TRUE setting flag to 1\n";
	  $bindingflag=1;
	  %selectedbinding=%bindingentry;
	  $arrayref=$j;	
	 }
	else {
		#print "flag is $bindingflag\n";
		};

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

sub smallest
{
my @array=@_;
@array=sort {$a<=>$b} @array;
return $array[0];
}

sub scopsearchhash
{
my @search=@{$_[0]};
my $cnt=@search;
my $lastid;
my $newid;
my %outhash;

for (my $i=0; $i<$cnt; $i++)
{
my $line=$search[$i];
$newid=substr($line,1,4);
#print "$newid=$lastid";
if ($newid ne $lastid)
{
#print" false\n";
#print "New id $newid index:$i\n";
$outhash{$newid}=$i;
$lastid=$newid;

}
elsif ($newid eq $lastid)
 {
# print " true\n";
}
}
return %outhash;


}


1;