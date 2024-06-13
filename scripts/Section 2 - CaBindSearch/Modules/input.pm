sub parsexml   #parses the data from spasm
{
my $openfile=$_[0];
print "location of file :$openfile\n";
open(FILE,"<$openfile");
my @readfile = <FILE>;
close FILE;

#print "first line :$readfile[0]\n";

 if ($readfile[0] eq "<XML>\n")
 {
  print "xml file found\n";
  if($readfile[1] eq "<smh type='hits'>\n")
  {
    print "xml type = hits\n";
   @outarray=parsehitxml(\@readfile);
   return @outarray;
  }
 elsif($readfile[1] eq "<smh type='bhits'>\n")
  {
    print "xml type = hits\n";
   @outarray=parsehitxml(\@readfile);
    return @outarray;
  }
elsif($readfile[1] eq "<smh type='mhits'>\n")
  {
    print "xml type = hits\n";
   @outarray=parsehitxml(\@readfile);
    return @outarray;
  }
   elsif($readfile[1] eq "<smh type='families'>\n")
  {
    print "xml type = families\n";
   @outarray=parsefamilyxml(\@readfile,"f");
   return @outarray;
  }
    elsif($readfile[1] eq "<smh type='rms'>\n")
  {
    print "xml type = rms\n";
   %outhash=parsermsxml(\@readfile);
return %outhash;
  }
elsif($readfile[1] eq "<smh type='data'>\n")
  {
    print "xml type = data\n";
   @outarray=parsefamilyxml(\@readfile,"d");
    return @outarray;
  }
elsif($readfile[1] eq "<smh type='input'>\n")
  {
    print "xml type = input\n";
   @outarray=parsefamilyxml(\@readfile,"i");
    return @outarray;
  }
 }
 else
 {
   print "no xml file found please check format and try again\n";
 }
 
}

sub parsehitxml
{
my @readfile=@{$_[0]};

my @outarray;

my $inhit='0';

my @res;

my $runcount;

print "parse hit xml initialised\n";

$linecount=@readfile;
#print "line count: $linecount\n";

for (my$i=3;$linecount>$i;$i++)
{

$line=$readfile[$i];
 
 if (!$inhit)
  {
  #print "$i not in hit\n";
    if ($line =~m/<hit id\=\w*?:.:\d{1,4}-\d{1,4}.*>\n/)
    #if ($line =~m/<hit id\=\w{4}\d{1,4}-\d{1,4}>\n/)
    {
    $inhit=1;
    $i--;
   #print "hit found:\n";
    }
   elsif ($line =~m/<runcount>(.*)<\runcount>\n/)
    #if ($line =~m/<hit id\=\w{4}\d{1,4}-\d{1,4}>\n/)
    {
    $runcount=$1;
    }
  }
 else
  {
  my @res;
  my %ligand;
  my @metals;
      while($inhit)
        {    

  
        $line=$readfile[$i];
        #print "$i is in hit\n";
	if ($line =~m/<hit id\=(\w*?:.:\d{1,4}-\d{1,4})>\n/)
        #if ($line =~m/<hit id\=(\w{4}\d{1,4}-\d{1,4})>\n/)
          {
          $hitid=$1;
          $inhit=1;
          $hitcount++;
          print "hitid: $hitid\n";
          }
	if ($line =~m/<hit id\=(\w*?:.:\d{1,4}-\d{1,4}) run=(.*)>\n/)
        #if ($line =~m/<hit id\=(\w{4}\d{1,4}-\d{1,4})>\n/)
          {
          $hitid=$1;
	  $run=$2;
          $inhit=1;
          $hitcount++;
          print "hitid: $hitid\n";
          }
        elsif ($line =~m/<chain>(.*)<\\chain>\n/)
          {
          $chain=$1;
         # print "Chain found:$chain\n"; 
          }
        elsif ($line =~m/<binds>(.*)<\\binds>\n/)
          {
          $binds=$1;
	 my %met;
	$met{met}=$binds;
	push (@metals,\%met);
          #print "Binding found:$binds\n"; 
          }
        elsif ($line =~m/<scop>(.*)<\\scop>\n/)
          {
          $scop=$1;
          #print "Scop found:$scop\n";
          }
        elsif ($line =~m/<res>(\d*)<\\res>\n/)
          	{
		if($inligand)
			{
			my $lires=$1;
	  		$ligand{res}=$lires;
			}
		else
			{
          		my$res=$1;
          		push (@res,$res);
          		print "Res found:$res\n";
          		#print "resarray:@res\n";
			}
         }
	elsif ($line =~m/<motif>(.*)<\\motif>\n/)
		 {
          $motif=$1;
         
          }
	elsif ($line =~m/<motif>(.*)<\\motif>\n/)
		 {
          $motif=$1;
         }
		
	
        elsif ($line =~m/<\\hit>/)
          {
         # print "end of hit found\n\n";
	  $hitid=~m/(\w*?):.:\d{1,4}-\d{1,4}/;
          my $pdbid=$1;
	  
          my%hit=("pdbid",$pdbid,"file",$hitid,"chain",$chain,"scop",$scop, "pdbin", $run,"motif",$motif,"bvalue",$bvalue);
          $hit{resnum}=\@res;
	  $hit{ligand}=\%ligand;
	  $hit{metal}=\@metals;
          push(@outarray,\%hit);
          $inhit=0;
          }

	elsif ($line =~m/<ligand>(.*)<\\motif>\n/)
		{
          	$inligand=1;
		%ligand=null;
         	}

	elsif ($line =~m/<resid>(.*)<\\resid>\n/&& $inligand)
          {
          my $liresid=$1;
	  $ligand{resid}=$liresid;
          }
        elsif ($line=~m/<bindvalue>(.*)<\\bindvalue>\n/)
	  {
	  $bvalue=$1;
	  }
	elsif ($line =~m/<\\ligand>\n/&&$inligand)
	{
          $inligand=0;
        }
        $i++;
      }
  }


}
my $count = @outarray;
print "$count hits found...\n";
if ($runcount!=0)
	{
	my %resultshash=("runcount",$runcount,"hitarray",\@outarray);
	return %resultshash;
	}

return @outarray;


}

sub parsefamilyxml
{
              print "parse family xml initialised\n";

my@readfile=@{$_[0]};
my $typeflag=$_[1];
my $tothits;

my@outarray;

my$inscop='0';

my $linecount=@readfile;
#print "line count: $linecount\n";

for (my$i=3;$linecount>$i;$i++)
{
$line=$readfile[$i];

 if (!$inscop)
  {
 # print "$i not in scop\n";
    #if ($line =~m/<hit id\=(\w{4}:.:\d{1,4}-\d{1,4})>\n/)
    if ($line =~m/<scopclass id\=\w\.\d{1,3}\.\d{1,3}\.\d{1,3}>\n|<scopclass id\=not found>\n/&& $typeflag eq "f")
    {
    $inscop=1;
    $i--;
  #  print "scopclass found:\n";
    }
    elsif ($line =~m/<scopclass id\=\w\.\d{1,3}\.\d{1,3}>\n|<scopclass id\=not found>\n/&& $typeflag eq "d")
    {
    $inscop=1;
    $i--;
 #   print "scopclass found:\n";
    }
    elsif (($line =~m/<scopclass id\=\w\.\d{1,3}\.\d{1,3}>\n|<scopclass id\=not found>\n/||$line =~m/<scopclass id\=\w\.\d{1,3}\.\d{1,3}\.\d{1,3}>\n|<scopclass id\=not found>\n/)&& $typeflag eq "i")
    {
    $inscop=1;
    $i--;
 #   print "scopclass found:\n";
    }
  }
 else
  {
  my @pdblist;
      while($inscop)
        {
	
        $line=$readfile[$i];
        #print "$i is in hit\n";
        #if ($line =~m/<scopclass id\=(\w\.\d{1,3}\.\d{1,3}\.\d{1,3})>\n|<scopclass id\=(not found)>\n/)
        if ($line =~m/<scopclass id\=(.*)>\n/)
          {
          $scopid=$1;
          $hitcount++;
 #         print "scopid: $scopid\n";
          }
          elsif ($line =~m/<id>(.*)<\\id>\n/)
          {
          my$pdb=$1;
	  print "*PDB*:$pdb\n";
	  $tothits++;
          push (@pdblist,$pdb);
          print "pdbid found:$pdb\n";
   #       print "resarray:@pdblist\n";
          }
        elsif ($line=~m/<bindvalue>(.*)<\\bindvalue>\n/)
	{
	$bvalue=$1;
	}
       elsif ($line =~m/<\\scopclass>/)
          {         
    #      print "end of hit found\n\n";
          my%scop=("scop",$scopid,"pdblist",\@pdblist,'bvalue',$bvalue);
          push(@outarray,\%scop);
          $inscop=0;
          }
        $i++;
      }
  }


}
my $familycnt=@outarray;
print "$familycnt families found...\n";
print"with $tothits total hits...\n";
return @outarray;  
             
           }


sub addxmlscopdata
{
 my @hitlist=@{$_[0]};
 my $filename=$_[1];

 my @scoplist=parsexml($filename);

 my $hitcnt=@hitlist;
 my $scopcnt=@scoplist;
 print "addscop data\n";
 for (my $x=0;$scopcnt>$x;$x++)
 {
   my %scop=%{$scoplist[$x]};
   my $id=$scop{scop};
   #print "scop id :$id\n";
   my @pdbarray=@{$scop{pdblist}};
   for $eachpdb(@pdbarray)
   {
     #print "pdb:$eachpdb\n" ;

     for (my $y=0;$hitcnt>$y;$y++)
     {
       my %eachhit=%{$hitlist[$y]};
	my $eachhitfile=$eachhit{file};	
      if ($eachpdb eq  $eachhitfile)
      {
      #print "$id added to hit $eachhit{file}\n";
      $eachhit{scop}=$id;
      $hitlist[$y]=\%eachhit;
      }
	else
	{
	#print "$eachpdb not equal to $eachhitfile\n";
	}
     }
   }
   
 }
return @hitlist;   
}

sub parsermsxml
{
 #print "parse rms xml initialised\n";
@file=@{$_[0]};
  my %rmshash;
	for my$line (@file)
		{
		if($line=~m/<RMS ID=.*\sID=.*>\s.*\s<\\RMS>/)
{
		$line=~m/<RMS ID=(.*)\sID=(.*)>\s(.*)\s<\\RMS>/;
		my $ID1=$1;
		my $ID2=$2;
		my $rms=$3;
		#print" line read : id:$ID1, ID:$ID2, RMS:$rms\n";
		
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
my @rmshash=%rmshash;
print "finished hash table: @rmshash\n";
  return %rmshash;
           }

sub addxmlrmsdata
{
 my @hitlist=@{$_[0]};
 my $filename=$_[1];

 my %rmslist=parsexml($filename);

my $hitlistcnt=@hitlist;

for (my $i=0; $i<$hitlistcnt; $i++)

{
my %hit=%{$hitlist[$i]};
my $hitid=$hit{file};
#print"hit = $hitid\n";
my%rmsdata=%{$rmslist{$hitid}};
my@rmsdata=%rmsdata;
#print"matches : @rmsdata\n";

$hit{match}=\%rmsdata;
$hitlist[$i]=\%hit;
}
return @hitlist;
}














1;