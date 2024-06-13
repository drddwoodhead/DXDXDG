sub iterativerun #handles the looping of spasm takes an array of hit ids as input
{
#print "\nITERATIVE RUN CALLED\n\n";

my $database=$_[1];
my @trylist=@{$_[0]};
my $sitedir=$_[3];
my $ftpdir=$_[2];
#$dbfile='spasm_25_jan07.lib';
$database=~m/(^.*\/)((.+)\.lib)/;
#print "1:$1\n2:$2\n3:$3\n";
my $sdb=$3;
my $spasmdatabase=$2;
my $spasmdir=$1;


my $rmscutoff=$_[4];
my $acmiss=$_[5];
my $scmiss=$_[6];
my $logdir=$_[7];
my $resultsdir=$_[8];
my $runident=$_[9];
my $runcount=$_[10];
my @runlist=@{$_[11]};
my $iteration=$_[12];
my $trylistcnt=@trylist;
#my $trylistcnt=2;
#my $runlistcnt=@runlist;
#print "sitedir*:$sitedir\n";
for (my $t=0;$t<$trylistcnt;$t++)
	{
	my $tryfileid=$trylist[$t];
	#print "tryfileid $tryfileid\n";
	$tryfileid=~m/(\w{4}):(\w):(\d{1,4})-(\d{1,4})/;
  	my $tryfile=$sitedir.$1."_".$2."_"."$3\-$4\.pdb";  
        $runflag=1;
	#print "for each entry in the try list (master hit list)\n";
	my $run=$t+1;
	print "\n\nrun $run of $trylistcnt\n";
	my $runlistcnt=@runlist;
	for ($r=0;$r<$runlistcnt;$r++)
		{
		$run=$runlist[$r];
		
		#print "check against hits already run\n";
		

		if ($tryfileid eq $run)
			{
			#print "$r of $runlistcnt try is equal to run $tryfileid = $run\n";
			
			$runflag=0;
			}
		else
			{
			 #print "$r of $runlistcnt try is notequal to run $tryfileid \/= $run\n";
			}
		}

if ($runflag)
			{  
			my @spasmresult=runspasm($tryfile,$database, $run, $rmscutoff, $logdir,$spasmdir,$rmscutoff,$acmiss,$scmiss);
			my $pdbfilename=$tryfile; $pdbfilename=~s/\/.*\///;
#  			print "Pdbfile:$pdbfilename\n";
			my $runcnt=@runlist;
			print "adding $tryfileid to tried list, $runcnt runs counted\n ";
   			push(@runlist,$tryfileid);
   			print "spasm finished, parsing results\n";
			#print "spasmresult: @spasmresult\n";
			if(@spasmresult)
			{

			@currenthitarray=parsespasm($j,\@spasmresult,$iteration);
			undef (@spasmresult);
			
			my $currentcnt=@currenthitarray;
			#print "$currentcnt hits\n";

			#print "getting pdb files\n";
   			#@ftp=`$ftpscript $macro`;

			$runlistcount=@runlist;
			@trylist=addresults(\@trylist,\@currenthitarray,$ftpdir,$sitedir,$resultsdir,$runident,$runlistcount);
			undef (@currenthitarray);
			}
			}
}
my %lists=(trylist,\@trylist,runlist,\@runlist);
return %lists;			
}

sub runspasm #runs spasm first checks to see if the run has been attempted before and is recorded in a log file takes a pdb file location as an input. flag can force the run to be redone.
{ 
  my $database=$_[1];
  #print "database:$database\n";
  my $rmscutoff=$_[3];
  my $run=$_[2];
  my $logdir=$_[4];
  my $spasmdir=$_[5];
  my $rmscutoff=$_[6];
  my $acmiss=$_[7];
  my $scmiss=$_[8];
 # print "*runspasm called...\n";
  my $pdbfile=$_[0];
  $pdbfile=~m/^.*\/(.+?).pdb/;
  my $pdbid=$1;
  #print "PDBFILE:$pdbfile\n";
  print "running spasm on $pdbid\n";
  my $filename="$logdir$pdbid\.log";
  #print "LOG:$filename\n";
  

 my $tempfile=$spasmdir."temp.pdb";
  open (PDBFILE,"<$pdbfile");
  my @pdbfile=<PDBFILE>;
  close PDBFILE;
  my $linecnt=@pdbfile;
	for (my $i=0;$i<$linecnt;$i++)
	{
	my $line=$pdbfile[$i];
	#print "$line";
	if ($line=~m/^HETATM.*/)
		{
		#print "deleting HETATM line: $line";	
		delete $pdbfile[$i];
		splice(@pdbfile, $i, 1);	
		}
	else {}
	}  

  open (TEMPFILE,">$tempfile");
	for my $line(@pdbfile)
	{print TEMPFILE $line;}
  close TEMPFILE;
  $database=~m/^.*\/(.+)\.lib/;
  my $db=$1;
  #print "DB:$db\n";
  #$spasmdatabase=$2;
  #$spasmdir=$1;

  my$inid=$pdbfile;
  $inid=~s/.*\///;


  my$matchflag=0;
  my@results=();


if ($forcespasm==0)
{
my $runcomplete=0;
#print "no force spasm!\n";
#print "Log:$filename\n";
 if (open(LOG,"<$filename")) 
{
	$runcomplete=0;

	my $filesize = -s "$filename";
   	#print "opening file:\n $filename\n...\n";
	if($filesize=>8000)
		{
		#print "file found to be over 8000b\n";
		$runcomplete=1;
		}
	else	{
		#print "file found to be LESS than 8000b\n";
		}
}
else {}

if ($runcomplete==1)
{
	open(LOG,"<$filename");
	@readlog=<LOG>;
	close LOG;

	$runcomplete=0;
	for my $line(@readlog)
	{
		$line=~ s/^ //;
		if ($line =~/.*Run again \? \(Y\)  Run again : \(N\).*/)
 			{
			#print "end of file found\n";
			$runcomplete=1;
			} 
		else {
			#print "$line\n";
			}
	}
}
else {}


if ($runcomplete==1)
	{
my $matchcount=0;
	my $linecnt=@readlog;
CHECK: for (my $i=0;$linecnt>$i;$i++)
		{
		$line=@readlog[$i];

		if ($line=~/^*MATCH #/)
      			{
			$matchcount++;
      			$matchflag=1;
      			while ($matchflag)
        			{
				$i++;
				$line=@readlog[$i];
               			if($line=~/^*Total BLOSUM-45/)
          				{
					$matchflag=0;
					}
				elsif($line=~/^*MATCH #|.*--------------------------------------------------.*|.*Run again \? \(Y\)  Run again : \(N\).*/)
					{
					print "$line\n";
					$runcomplete=0;
					$matchflag=0;
					#print "problem found in hit; file corrupt\n";
					last CHECK ;

					}
				else 
					{
					}
				}
			}

		else {}
		}
	#print "file read $matchcount hits found!\n";
	}
else {}

if ($runcomplete==1)
{
print "stored file loaded $filename\n";
#print "@readlog\n";
return @readlog;
}
else
{
#print "error logfile loaded but NOT complete\n Rerunning spasm\n";
}
}


elsif($forcespasm==1)
{print "force spasm\n";}
else{}


my $filecomplete=0;

if (open(INPUT,"<$pdbfile")) 
	{
	open(INPUT,"<$pdbfile");
	my @readin=<INPUT>;
	close INPUT;
	$readin=@readin;

	if (@readin)
 		{
		$filecomplete=1;
		} 
	else {}
	}
else{}

if(!$filecomplete)
	{
	$pdbid=~m/(.{4}):(.):(\d{1,4})-(\d{1,4})/;
	my $pdbcode=$1;
	my $chain=$2;
	my $res1=$3;
	my $res2=$3+2;
	my $res3=$4;
	my @arrayres=($res1,$res2,$res3);
	mkpdb($pdbcode, $ftpdir, \@arrayres, $chain);
	}
else{}


my $substitutions=findsubs($tempfile);

if ($substitutions)
	{
	
  	#print "substitutions used: \n*$substitutions*\n";
  	#print "$spasm $database $tempfile $inid $rmscutoff $acmiss $scmiss $substitutions\n";
  	$database=~s/\s*//g;
  	my @results=`$spasm $database $tempfile $inid $rmscutoff $acmiss $scmiss \"$substitutions\"`; #runs spasm script with location of database and pdb id as input
  	#print @results;
	#`rm $tempfile`;
  	print "SPASM finished\n";
  	#my$filename="spasm.log";
  	#print "saving spasm log file:\n $filename \n";
  	open(SPASMLOG,">$filename"); # open new log file
	for my$line (@results)
		{
		print SPASMLOG "$line";
		}


  	#print SPASMLOG @results;
  	close SPASMLOG;

  	return @results;
	}
else
	{
	return '0';
	} 

}

sub findsubs  #creates text string with appropiate substitutions for the residues to be fed back into spasm
{
 my $openpdb=$_[0];
 #print "open pdb : $openpdb\n";
 my @resarray;
 my $resnumber;
 my $returnres="";

  open(PDBFILE,$openpdb);
 @readpdb = <PDBFILE>;
 close PDBFILE;
 #print @readpdb;

 for $line(@readpdb)
  {
    if ($line =~ /END/)
    {}
    else {
  my $residue=substr($line,17,3);
  my $resfound=0;
  $resnumber=@resarray;
  for (my $i=0; $i<$resnumber;$i++)
    {

    $res=$resarray[$i]   ;
    #print "residue $residue\n";
    #print "$i record:$record:$pdb:pdb\n";
       if ($residue eq $res)
          {
            $resfound=1;
          }
       else{}
    }
  if (not $resfound)
   {
   #print "$residue being written and added resarray\n";
            push(@resarray, $residue);
   }
   }
   }
   @resarray=sort(@resarray);
   #print "res array: \n @resarray\n";
      $arraysize=@resarray;
      #print "arraysize $arraysize\n";

     for (my $j=0; $j<$arraysize;$j++)
    {
      my $subres=$resarray[$j];
      #print "$subres";
    #print "****res $subres **** \n";
	my $sub=substitute($subres);
if($sub)
{
    $returnres=$returnres.$sub;
 }
else
{return '0';}
   }
   chomp($returnres);
   #print "\$returnres:\n $returnres";
undef(@readpdb);

   return $returnres;
}

sub substitute # contains table of acceptable substitutions for use with findsubs
{
    my($residue) = @_;

    $residue = uc $residue;

    my(%substitutions) = (

    'ASN' => "ASN ASP SER\n",
    'ASP' => "ASP ASN SER\n",
    'SER' => "SER ASN ASP\n",
     );

    if(exists $substitutions{$residue}) {

        return $substitutions{$residue};
    }else{


            print STDERR "Bad residue \"$residue\"!!\n";
            return '0';
    }
}
 
sub parsespasm   #parses the data from spasm makes an array of hash tables each associated with a hit
{
    #print "*parsespasm called...\n";
    #print "@_\n";
    my @results=@{$_[1]};
    my $iteration=$_[2];
    #print @results;
    my $i=0 ;
    my @matcharray;
    
    $matchcount=0;
    $compflag=0;
    my $iter=$_[0];
    my $initialinfoflag=1;
    #print "iter:$iter\n";


  while ($results[$i])
    {
    $line=$results[$i];
    ++$i;

    if ($line =~ /^*Four-character ID for this run/)
    {

    $line=~m/^\s{0,2}Four-character ID for this run \? \(.{4}\)  Four-character ID for this run : \((.{4})\) /;
    $runid=$1;
    }
  
    elsif ($line =~ /^*==> HIT/)
      {
      $line=~m/^\s{0,2}==> HIT  : \((.{4})\)/ ;
      $line=$1;
      $line=lc($line);
      $hittable{pdbid}=$line;
      #$hittable{iter}=$iter;
      
      }
    elsif ($line =~ /^*Compound/)
      {
      $line=~s/Compound//g;$line=~s/\s//g;$line=~s/\(//;$line=~s/\)//;$line=~s/\://;
      $hittable{compound}=$line;
      $compflag=1;
      }
    elsif ($line =~ /^*Entry/)
      {
        if($compflag){$compflag=0;}
      }
    elsif ($line=~/^*MATCH #/)
      {
      $absolutehits++ ;
      $rms=substr($line,25,4);
#print "rms value recorded (25,4):$rms\n";
     
      my %match=%hittable;
      my @resarray;
      my @restype;
       my @mresarray;
     my @mrestype;
     
      $matchflag=1;
      while ($matchflag)
        {
        $line=$results[$i];

         ++$i;
        if($line=~/^*Total BLOSUM-45/)
          {
       		              
                my $count=@resarray-1;
		my$fileid=$match{pdbid}.":".$match{chain}.":".$resarray[0]."-".$resarray[$count];
		$match{file}=$fileid;
        #print "this match file id is :: $match{file}";
        $match{resnum}=\@resarray;
          $match{restype}=\@restype;
		$match{bvalue}="B+";
		$match{iteration}=$iteration;
	
	
	#my $mcount=@mresarray-1;
	#my$rmsid=$runid.":".$match{chain}.":".$mresarray[0]."-".$mresarray[$mcount];
	#my%rms;
	#$rms{$rmsid}=$rms;
	
	#print "RMS hash (parser):\n ::"."@{[ %rms ]}\n"."\n::\n\n";
	
	#$match{match}=\%rms;
	
          	push(@matcharray,\%match);
          	#print ("pdbid:", $matcharray[$matchcount]{pdbid}, "\n")  ; #to access the pdb id
          	#print ("residue matches:", "@{$matcharray[$matchcount]{resnum}}", "\n"); #to access array of residues
          	#print ("RMS value:", $matcharray[$matchcount]{rms}, "\n");
          	++$matchcount;
	
          $matchflag=0;
          }
        else
          {
	  $line=~m/^\s{0,2}\w{3}\s\w\s{0,3}\d{1,4}\s{0,4}<--->\s\w{3}\s(\w)\s{0,3}\d{1,4}/;
          my $chain=$1;
		if(!$chain||($chain eq " ")){$chain='-';}
          $match{chain}=$chain;
	#print "line $line\n chain $chain\n";

          #my$res= substr($line,23,4);
	  
	  $line=~m/^\s{0,4}.{11}<--->.{6}\s{0,3}(\d{1,4})/;
	  my $res=$1;
          $res=~s/\s//g;
	  $res=~s/\D//g;	
          #print ("residue number:$line\n");
          push(@resarray,$res);

           my$mres= substr($line,7,3);
          $mres=~s/\s//g;
          #print ("residue number:$line\n");
          push(@mresarray,$mres);

           my$rest= substr($line,18,3);
          $rest=~s/\s//g;
          #print ("residue number:$line\n");
          push(@restype,$rest);

          # my$mrest= substr($line,1,3);
          #$res=~s/\s//g;
          #print ("residue number:$line\n");
          push(@mrestype,$mrest);
          }
        }


      my $totres=@mresarray;
      #print "runid: $runid\n";

      $rmsid=$runid.":".$chain.":".$mresarray[0]."\-".$mresarray[$totres-1];
      $rmsid=~s/\s//g;
	
      #$rmsid
      #print "this RMS has been recorded $rmsid\n";
      #print "rms:$rms\n";
      $match{$rmsid}=$rms; #*****************************
      $match{pdbin}=$runid;#*******************************

      }
    else
      {
      if($compflag)# if the compound description is on multiple lines this will capture it
        {
        $line=~s/\)//;
        $line=~s/\n//;
        $hittable{compound}= $hittable{compound}.$line;
        }
      }
    }
undef(@results);
my $count=@matcharray;
print "$count hits found in spasm output\n";
#print"removing short hits (other than $resexpect residues...\n";
@matcharray=removeshort(\@matcharray,$resexpect);
$count=@matcharray;
print "$count hits\n";

return @matcharray;
}

sub addresults #adds the results from the current iteration of spasm search to the hits file, makeing
               #sure that only one copy of each hit is recorded. takes an array of already present hit ids and the new array of hits as inputs.
{
 #print "*addresults called...\n";
 my @oldhitlist=@{$_[0]};
 my @newhits=@{$_[1]};
 my $ftpdir=$_[2];
 my $sitedir=$_[3];
 my $resultsdir=$_[4];
 my $runident=$_[5];
 my $runcount=$_[6];

#print "sitedir**:$sitedir\n";

my @newhitlist;#contains all the new hits that have not already been found.

 my $newhitcnt=@newhits;
 my $oldhitcnt=@oldhitlist;
 my $found=0;
 #my @outhits=@oldhits;
 #print "@oldhits,@newhits";
 #print "NEWHITSIN: @newhits\n";
 #print "OLDHITSIN: @oldhitlist\n";

 for ($y=0; $y<$newhitcnt; $y++)
        {
        $found=0;
        my%new=%{$newhits[$y]};

        for ($x=0; $x<$oldhitcnt; $x++)
            {
            
            if ($new{file} eq $oldhitlist[$x])
             	 {
		#print "$new{file} == $oldhitlist[$x]\n";
		$found=1;
	     	 }
	    else
		 {
		#print "$new{file} != $oldhitlist[$x]\n";
                #print "not equal newfile:$new{file},oldfile:$old{file}\n";
           	 }
	          }

        if (!$found)
	          {
            #print "############################\nnew hit found push  $new{pdbid} @{$new{resnum}} \n";
             push (@newhitlist, \%new) ;
		#print "NEWHITLIST:@newhitlist\n";
	     
            }
    	  }
#if ($e!=$iterations-1)
#{
$hitcount=@newhitlist;
print "$hitcount hits\n";

my $newhitcnt=@newhitlist;
print "making pdb files...\n";
for (my $i=0; $i<$newhitcnt; $i++)
        			{
       				my %each=%{$newhitlist[$i]};
       				my $ref=\@{$each{resnum}} ;
       				#print "making pdb ".$each{pdbid}." residues "."@{$ref}"."\n";
				mkpdb ($each{pdbid}, $ftpdir, $ref, $each{chain},$sitedir);
       				}

print "adding motif data...\n";
@newhitlist=addmotifdata(\@newhitlist,$ftpdir);
##$hitcount=@newhitlist;
$hitcount=@newhitlist;
print "$hitcount hits\n";

print "sort by motif data...\n";
@newhitlist=sortmotifdata(\@newhitlist);	
$hitcount=@newhitlist;
print "$hitcount hits\n";

print "finding bound metals...\n";
@newhitlist=addmetalbinding(\@newhitlist,$ftpdir,$sitedir);
$hitcount=@newhitlist;
print "$hitcount hits\n";

print "removing hits that don't bind anything...\n";		
@newhitlist=findbound(@newhitlist);
$hitcount=@newhitlist;
print "$hitcount hits\n";
#print "adding extra ligand data and removing those with no extra ligand...\n";
#@newhitlist=addliganddata(\@newhitlist);
$hitcount=@hitlist;
$hitcount=@hitlist;
print "$hitcount hits\n";
#print "removing hits that don't bind anything...\n";		
#@newhitlist=findbound(@newhitlist);
#$hitcount=@newhitlist;
#print "$hitcount hits\n";
#print "adding extra ligand data and removing those with no extra ligand...\n";
#@currenthitarray=addliganddata(\@currenthitarray);
#$hitcount=@currenthitarray;
#$hitcount=@newhitlist;
#print "$hitcount hits\n";

print "saving to file...\n";
savehits(\@newhitlist,"a",$resultsdir,$runident,$runcount);
$newhitcount=@newhitlist;
print "$newhitcount new hits found...\n";

for(my$i=0;$i<$newhitcount;$i++)
{
%new=%{$newhitlist[$i]};
push (@oldhitlist, $new{file}) ;
}

undef(@newhitlist);
#print "OLDHITS: @oldhitlist\n";
return @oldhitlist;
}

sub initialisearray
{
my $openpdb=$_[0];
my $dir=$_[1];

$openpdb=~m/^.*?\/(.*)/;
my $id=$1;
$id=~s/_/:/g;
my @outarray;
my %first;


#my $openpdb ="$dir$id.pdb";
#print "Initialising Try List...\nLocatoin of initial pdb :$openpdb\n";
 open(PDBFILE, $openpdb);
my  @readpdb = <PDBFILE>;
 close PDBFILE;

 #print "@readpdb\n";
 
my@resarray;
my@restypearray;

for my $line(@readpdb)
{

my $res=substr($line,22,4);
 $res=~ s/\s//g;
my $restype=substr($line,17,3);
my $flag=0;
  if ($resarray[0]eq ""||!@resarray)
    {     
	#  print "res: $restype $res \n";
    push(@restypearray,$restype);
	  push(@resarray,$res);
    }

  else
    {
	  for $resfound(@resarray)
		  {
		  if($res eq $resfound)
		      {$flag=1;}
      else
      {}
        
	    }
      if (!$flag&&$res)
        {
	      #  print "res: $restype $res \n";
          push(@restypearray,$restype);
	        push(@resarray,$res);
	        }	
    }
}
undef(@readpdb);

$id=~m/((.{4}):(.):\d{1,4}-\d{1,4})/;
my $chain=$3;
if(!$chain){$chain="-";}

my %rms=($id,'0');

$first{file}=$1;
$first{filepath}=$id;
$first{chain}=$chain;
$first{pdbid}=$2;
$first{resnum}=\@resarray;
$first{restype}=\@restypearray;
$first{match}=\%rms;
$first{iteration}="0";
#print "firstid : $first{pdbid}\n";
#print "seedpdb : $first{file}\n";
#print "chain   : $first{chain}\n";
#print "res : ".@{$first{resnum}}."\n";


push (@outarray, \%first);
$trylistcnt++;
return @outarray;
}

sub removeshort #removes hits with anything other than the expected number of residues from an array of hits
{
my @hitlist=@{$_[0]};
my $resexpect=$_[1];
my @outarray;
my$hitlistcnt=@hitlist;

for (my $i=0;$hitlistcnt>$i;$i++)
{
my %hit=%{$hitlist[$i]};
#print "metal bound:$hit{metal}\n";

my @resarray=@{$hit{resnum}};
my $rescheck=@resarray;

#print ":@resarray";

if($rescheck==$resexpect)
{
#print "metal bound :true\n";
push(@outarray,\%hit);
}



}
return @outarray;
}


1;