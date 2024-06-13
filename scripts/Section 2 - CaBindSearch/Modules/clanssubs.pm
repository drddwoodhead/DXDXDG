
sub mkseqfile  #makes a file containing the sequences in fasta format
{
my@pdblist=@{$_[0]};
my $pdbdir="/home/dwoodhead/wrkdir/ftpdir/";

#$filename ="seqences\_$seedpdb\_$dbfile";
my $filename=$resultsdir.$runident."-c\.fas";
open(SEQFILE,">$filename"); # open new sequence file
$pdblist=@pdblist;

for (my $i=0; $i<$pdblist; $i++)
  {
  my $seqdone=0;
  my $pdb=lc($pdblist[$i]{pdbid});
  my %selectpdb=%{$pdblist[$i]};
  my $pdb=$selectpdb{scop}."_".$selectpdb{pdbid};
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
  $pdbid=substr($pdb,-4,4);
  my$openpdb ="$pdbdir$pdbid.pdb";
  open(PDBFILE, $openpdb);
  @readpdb = <PDBFILE>;
  
  # Parse the record types of the PDB file
  my %recordtypes = parsePDBrecordtypes(@readpdb);
 close PDBFILE;
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

 sub brokenmkmotifseqfile  #makes a file containing the motif sequences in fasta format
{
my@pdblist=@{$_[0]};
my$range=$_[1];
my $pdbdir="/home/dwoodhead/wrkdir/ftpdir/";

my $filename=$resultsdir.$runident."motif\.fas";
open(SEQFILE,">$filename"); # open new sequence file
$pdblist=@pdblist;

for ($i=0; $i<$pdblist; $i++)
  {
  my $seqdone=0;
  my %selectedhit=%{$pdblist[$i]};
  my $pdb=lc($selectedhit{pdbid});
  my $pdbid= $selectedhit{scop}."_".$selectedhit{file} ;
 # print "$pdb\n";
  my @selectedres=@{$selectedhit{resnum}};
  my $restot= @selectedres;
  my $res1= $selectedres[0] ;
  my $res2= $selectedres[$restot-1] ;
  my $seqnumber=@seqrecord;
  my $chain=$selectedhit{chain};
  #print "chain:: $chain\n";
  #print "pdb:$pdb\n";
  my$openpdb ="$pdbdir$pdb.pdb";
  open(PDBFILE, $openpdb);
  @readpdb = <PDBFILE>;
  close PDBFILE;

  #print " Reabpdb : \n @readpdb\n";
  
  # Parse the record types of the PDB file
  my %recordtypes = parsePDBrecordtypes(@readpdb);

  # Extract the amino acid sequences of chosen chain in the protein
  my $selectedchain = extractSEQRESCHAIN( $recordtypes{'SEQRES'}, $chain );

  #print"Chains::: $selectedchain\n"   ;

 # Translate the 3-character codes to 1-character codes, and print
  print SEQFILE ">$pdbid\n";
  #foreach my $selectedchain (@chains)
    #{
    #print ">>>>$pdb\n";
    #print "****chain $chain **** \n\n";
    #print "$chain\n";
     $oneaacode=iub3to1($selectedchain);



      $startresidue=getstartresidue($openpdb,$chain) ;

     if($range)  
                   {
      my $lengthadjust=0;
      my $substrstart=$res1-$startresidue-$range+1;
		if ($substrstart<0)
			{
			$lengthadjust=$substrstart;
			$substrstart=0;
			}
    #   print "res1:$res1 \- startres: $startresidue \- Range:$range \= substart:$substrstart\n";
      my $substrlength=(($res2-$res1)+1)+(2*$range)-$lengthajdust;
     # print "\(res2:$res2 \- res1 $res1\) \+ \(2xrange:$range\) \= $substrlength\n";
      $oneaacode=substr($oneaacode,$substrstart,$substrlength);
	 print SEQFILE  "pdb start res:$startresidue, firstresidue:$res1, substring start:$substrstart\n";
                    }
     print SEQFILE  $oneaacode;
   # }
   print SEQFILE "\n\n";
  }
 close SEQFILE;

 return;
 }

sub mkmotifseqfile  #makes a file containing the motif sequences in fasta format
{
my@pdblist=@{$_[0]};
my$range=$_[1];


my $filename=$resultsdir.$runident."-r-$range-\.fas";
open(SEQFILE,">$filename"); # open new sequence file
$pdblist=@pdblist;

for ($i=0; $i<$pdblist; $i++)
  {
  my %selectedhit=%{$pdblist[$i]};
  
  my $pdbid= $selectedhit{scop}."_".$selectedhit{file} ;
  print SEQFILE ">$pdbid\n";  
  print ">$pdbid\n";

  my $oneaacode=getseq(\%selectedhit,$range);
  print SEQFILE  $oneaacode;
  print  $oneaacode
  # }
  print SEQFILE "\n\n";
  }
 close SEQFILE;

 return;
 }

sub rmsmatrix #produces a matrix file suitable to give clans as input
{
my @hitlist=@{$_[0]};


#$filename ="seqences\_$seedpdb\_$dbfile";
my$filename=$resultsdir."rmsmatrix-".$runident."\.clans";
open(MATRIXFILE,">$filename"); # open new sequence file
my$hitlist=@hitlist;

print MATRIXFILE "sequences=".$hitlist."\n";
print MATRIXFILE "<seqs>\n" ;
for (my $x=0; $x<$hitlist; $x++)
{
	my%hit=%{$hitlist[$x]};


print MATRIXFILE ">".$hit{scop}."_".$hit{file}."\n";
}
print MATRIXFILE "</seqs>\n";
print MATRIXFILE "<mtx>\n";

for (my $i=0; $i<$hitlist; $i++)
        {
	my%hitrow=%{$hitlist[$i]};
	#$hitrow=@hitlist[$i];
	#%rowrms=%{$hitrow{match}};
        #my $restotr=@{$hitrow{resnum}};
        $rowhitinfo=$hitrow{file};
        for (my $j=0; $j<$hitlist; $j++)
         	{
		#@hitcolumn=@hitlist[$j];
		my%hitcolumn=%{$hitlist[$j]};
    my @colresarray=@{$hitcolumn{resnum}};
    my $restotc=@colresarray;
		my $colhitinfo=$hitcolumn{file};
                #print "colhitinfo: $colhitinfo\n";
		#print "rowhitinfo: $rowhitinfo\n\n";
		my %rmstable=%{$hitrow{match}};

		if ($rowhitinfo eq $colhitinfo)   #hit colomn is the same as the row the difference is 0
                 {print MATRIXFILE "0 "}

		elsif($rmstable{$colhitinfo})    #the row has an entry under this colomns hit print it out
                 {print MATRIXFILE $rmstable{$colhitinfo}." ";}

		else              # if no hit is found and they are not the same hit put 5 as arbitary
                 {print MATRIXFILE "5 ";}
                 }
        print MATRIXFILE "\n";
        }
print MATRIXFILE "<\mtx>";

}

sub rmsmatrix2 #produces a matrix file suitable to give clans as input
{
my @hitlist=@{$_[0]};
my %rmshash=%{$_[1]};

#$filename ="seqences\_$seedpdb\_$dbfile";
my$filename=$resultsdir."rmsmatrix-".$runident.$rmsatomtypes."\.clans";
open(MATRIXFILE,">$filename"); # open new sequence file
my$hitlist=@hitlist;

print MATRIXFILE "sequences=".$hitlist."\n";
print MATRIXFILE "<seqs>\n" ;
for (my $x=0; $x<$hitlist; $x++)
{
my%hit=%{$hitlist[$x]};
print MATRIXFILE ">".$hit{scop}."_".$hit{file}."\n";
}
print MATRIXFILE "</seqs>\n";
print MATRIXFILE "<mtx>\n";

for (my $i=0; $i<$hitlist; $i++)
        {
	my%hitrow=%{$hitlist[$i]};
	#$hitrow=@hitlist[$i];
	#%rowrms=%{$hitrow{match}};
        #my $restotr=@{$hitrow{resnum}};
        $rowhitinfo=$hitrow{file};
        for (my $j=0; $j<$hitlist; $j++)
         	{
		#@hitcolumn=@hitlist[$j];
		my%hitcolumn=%{$hitlist[$j]};
    my @colresarray=@{$hitcolumn{resnum}};
    my $restotc=@colresarray;
		my $colhitinfo=$hitcolumn{file};
                #print "colhitinfo: $colhitinfo\n";
		#print "rowhitinfo: $rowhitinfo\n\n";
		
		my %rmstable=%{$rmshash{$rowhitinfo}};

		if ($rowhitinfo eq $colhitinfo)   #hit colomn is the same as the row the difference is 0
                 {print MATRIXFILE "0 "}

		elsif($rmstable{$colhitinfo})    #the row has an entry under this colomns hit print it out
                 {
		my $num=$rmstable{$colhitinfo};
		my $lognum=log($num);
		print MATRIXFILE $lognum." ";
		}

		else              # if no hit is found and they are not the same hit put 5 as arbitary
                 {print MATRIXFILE "0 ";}
                 }
        print MATRIXFILE "\n";
        }
print MATRIXFILE "</mtx>";

}

sub recipmatrix #produces a matrix file suitable to give clans as input
{
my@hitlist=@_;


#$filename ="seqences\_$seedpdb\_$dbfile";
my$filename="/home/dwoodhead/wrkdir/results/"."matrix\-".$seedpdb."\-".$dbfile."\-".$rmscutoff."\-".$iterations."\.fas";
open(MATRIXFILE,">$filename"); # open new sequence file
my$hitlist=@hitlist;

print MATRIXFILE "sequences=".$hitlist."\n";
print MATRIXFILE "<seqs>\n" ;
for ($x=0; $x<$hitlist; $x++)
{
	my%hit=%{$hitlist[$x]};


print MATRIXFILE ">".$hit{scop}."_".$hit{file}."\n";
}
print MATRIXFILE "</seqs>\n";
print MATRIXFILE "<mtx>\n";

for (my $i=0; $i<$hitlist; $i++)
        {
	my%hitrow=%{$hitlist[$i]};
	#$hitrow=@hitlist[$i];
	#%rowrms=%{$hitrow{match}};
        #my $restotr=@{$hitrow{resnum}};
        $rowhitinfo=$hitrow{file};
        for (my $j=0; $j<$hitlist; $j++)
         	{
		#@hitcolumn=@hitlist[$j];
		my%hitcolumn=%{$hitlist[$j]};
    my @colresarray=@{$hitcolumn{resnum}};
    my $restotc=@colresarray;
		my $colhitinfo=$hitcolumn{file};
                #print "colhitinfo: $colhitinfo\n";
		#print "rowhitinfo: $rowhitinfo\n\n";
		my %rmstable=%{$hitrow{match}};

		if ($rowhitinfo eq $colhitinfo)   #hit colomn is the same as the row the difference is 0
                 {print MATRIXFILE "20 "}

		elsif($rmstable{$colhitinfo})    #the row has an entry under this colomns hit print it out

                 {
		if($rmstable{$colhitinfo}==0)
		{
		print MATRIXFILE "20 "
		}
		else
		{
		my $reciprical=1/$rmstable{$colhitinfo};
		$reciprical=round($reciprical);
		print MATRIXFILE $reciprical." ";
		}
		}

		else              # if no hit is found and they are not the same hit put 5 as arbitary

                 {print MATRIXFILE "1 ";}
                 }
        print MATRIXFILE "\n";
        }
print MATRIXFILE "<\mtx>";

}
########################################################################################################################

sub extractSEQRES {

    use strict;
    use warnings;

    my($seqres) = @_;

    my $lastchain = '';
    my $sequence = '';
    my @results = (  );
    # make array of lines

    my @record = split ( /\n/, $seqres);

    foreach my $line (@record) {
        # Chain is in column 12, residues start in column 20
        my ($thischain) = substr($line, 11, 1);
        my($residues)  = substr($line, 19, 52); # add space at end

        # Check if a new chain, or continuation of previous chain
        if("$lastchain" eq "") {
            $sequence = $residues;
        }elsif("$thischain" eq "$lastchain") {
            $sequence .= $residues;

        # Finish gathering previous chain (unless first record)
        }elsif ( $sequence ) {
            push(@results, $sequence);
            $sequence = $residues;
        }
        $lastchain = $thischain;
    }

    # save last chain
    push(@results, $sequence);

    return @results;
}

sub extractSEQRESCHAIN {

    #use strict;
    #use warnings;

    my $seqres = $_[0];
    my $chain= $_[1];
    #print "seqreschain:$chain\n";
    my $lastchain = '';
    my $sequence = '';
    my @results = (  );
    # make array of lines

    my @record = split ( /\n/, $seqres);

    foreach my $line (@record)
    {
        # Chain is in column 12, residues start in column 20
        my ($thischain2) = substr($line, 11, 1);
        my($residues)  = substr($line, 19, 52); # add space at end
        #print "seqresthischain:::$thischain\n";
   
        if("$lastchain" eq "")#first chain found
        	{ $sequence = $residues;}
        elsif("$thischain2" eq "$lastchain")
        	{ $sequence .= $residues;}# Finish gathering previous chain (unless first record)
	elsif ( $sequence ) #new chain started
           	{
                if($lastchain eq $chain||$chain eq "-")
                	{
                	#print "$thischain equal $chain\n";
                	return $sequence;
                	}
		
             $sequence = $residues;
            	}
        $lastchain = $thischain2;
    
      }


      
  

    
}

sub parsePDBrecordtypes {

    my @file = @_;

    use strict;
    use warnings;

    my %recordtypes = (  );

    foreach my $line (@file) {

        # Get the record type name which begins at the
        # start of the line and ends at the first space

        # The pattern (\S+) is returned and saved in $recordtype
        my($recordtype) = ($line =~ /^(\S+)/);

        # .= fails if a key is undefined, so we have to
        # test for definition and use either .= or = depending
        if(defined $recordtypes{$recordtype} ) {
            $recordtypes{$recordtype} .= $line;
        }else{
            $recordtypes{$recordtype} = $line;
        }
    }

    return %recordtypes;
}

sub iub3to1 {

    my($input) = @_;

    my %three2one = (

    'ALA'=>'A', 'VAL'=>'V', 'PHE'=>'F', 'PRO'=>'P', 'MET'=>'M',
    'ILE'=>'I', 'LEU'=>'L', 'ASP'=>'D', 'GLU'=>'E', 'LYS'=>'K',
    'ARG'=>'R', 'SER'=>'S', 'THR'=>'T', 'TYR'=>'Y', 'HIS'=>'H',
    'CYS'=>'C', 'ASN'=>'N', 'GLN'=>'Q', 'TRP'=>'W', 'GLY'=>'G',
    '2AS'=>'D', '3AH'=>'H', '5HP'=>'E', 'ACL'=>'R', 'AIB'=>'A',
    'ALM'=>'A', 'ALO'=>'T', 'ALY'=>'K', 'ARM'=>'R', 'ASA'=>'D',
    'ASB'=>'D', 'ASK'=>'D', 'ASL'=>'D', 'ASQ'=>'D', 'AYA'=>'A',
    'BCS'=>'C', 'BHD'=>'D', 'BMT'=>'T', 'BNN'=>'A', 'BUC'=>'C',
    'BUG'=>'L', 'C5C'=>'C', 'C6C'=>'C', 'CCS'=>'C', 'CEA'=>'C',
    'CHG'=>'A', 'CLE'=>'L', 'CME'=>'C', 'CSD'=>'A', 'CSO'=>'C',
    'CSP'=>'C', 'CSS'=>'C', 'CSW'=>'C', 'CXM'=>'M', 'CY1'=>'C',
    'CY3'=>'C', 'CYG'=>'C', 'CYM'=>'C', 'CYQ'=>'C', 'DAH'=>'F',
    'DAL'=>'A', 'DAR'=>'R', 'DAS'=>'D', 'DCY'=>'C', 'DGL'=>'E',
    'DGN'=>'Q', 'DHA'=>'A', 'DHI'=>'H', 'DIL'=>'I', 'DIV'=>'V',
    'DLE'=>'L', 'DLY'=>'K', 'DNP'=>'A', 'DPN'=>'F', 'DPR'=>'P',
    'DSN'=>'S', 'DSP'=>'D', 'DTH'=>'T', 'DTR'=>'W', 'DTY'=>'Y',
    'DVA'=>'V', 'EFC'=>'C', 'FLA'=>'A', 'FME'=>'M', 'GGL'=>'E',
    'GLZ'=>'G', 'GMA'=>'E', 'GSC'=>'G', 'HAC'=>'A', 'HAR'=>'R',
    'HIC'=>'H', 'HIP'=>'H', 'HMR'=>'R', 'HPQ'=>'F', 'HTR'=>'W',
    'HYP'=>'P', 'IIL'=>'I', 'IYR'=>'Y', 'KCX'=>'K', 'LLP'=>'K',
    'LLY'=>'K', 'LTR'=>'W', 'LYM'=>'K', 'LYZ'=>'K', 'MAA'=>'A',
    'MEN'=>'N', 'MHS'=>'H', 'MIS'=>'S', 'MLE'=>'L', 'MPQ'=>'G',
    'MSA'=>'G', 'MSE'=>'M', 'MVA'=>'V', 'NEM'=>'H', 'NEP'=>'H',
    'NLE'=>'L', 'NLN'=>'L', 'NLP'=>'L', 'NMC'=>'G', 'OAS'=>'S',
    'OCS'=>'C', 'OMT'=>'M', 'PAQ'=>'Y', 'PCA'=>'E', 'PEC'=>'C',
    'PHI'=>'F', 'PHL'=>'F', 'PR3'=>'C', 'PRR'=>'A', 'PTR'=>'Y',
    'SAC'=>'S', 'SAR'=>'G', 'SCH'=>'C', 'SCS'=>'C', 'SCY'=>'C',
    'SEL'=>'S', 'SEP'=>'S', 'SET'=>'S', 'SHC'=>'C', 'SHR'=>'K',
    'SOC'=>'C', 'STY'=>'Y', 'SVA'=>'S', 'TIH'=>'A', 'TPL'=>'W',
    'TPO'=>'T', 'TPQ'=>'A', 'TRG'=>'K', 'TRO'=>'W', 'TYB'=>'Y',
    'TYQ'=>'Y', 'TYS'=>'Y', 'TYY'=>'Y', 'AGM'=>'R', 'GL3'=>'G',
    'SMC'=>'C', 'ASX'=>'B', 'CGU'=>'E', 'CSX'=>'C', 'GLX'=>'Z',
    'UNK'=>'X',
    'A'   => 'X',
    'C'   => 'X',
    'G'   => 'X',
    'T'   => 'X',
 #FGL , CSE , PVL, ACE, YCM, SNN  could not be found so placed as UNK
    'FGL'=>'X',
    'CSE'=>'X',
    'PVL'=>'X',
    'ACE'=>'X',
    'YCM'=>'X',
    'SNN'=>'X',


    );

    # clean up the input
    $input =~ s/\n/ /g;

    my $seq = '';

    # This use of split separates on any contiguous whitespace
    my @code3 = split(' ', $input);

    foreach my $code (@code3) {
        # A little error checking
        if(not defined $three2one{$code}) 
	{
            #print "Code $code not defined\n";
            #next;
	    $seq='X';
        }
	else{
        $seq .= $three2one{$code};
    		}
    }
    return $seq;
}

sub getstartresidue
{
my $openpdb=$_[0];
my $chain=$_[1];

#print "$openpdb\n";
open(PDBFILE, $openpdb);
my @readpdb = <PDBFILE>;
  
if ($chain && $chain ne "-" && $chain ne " ")
  {
 # print "chain info found $chain\n";
  my $cnt=@readpdb;
  for (my $i=0; $i<$cnt; $i++)
    {
    my $line=$readpdb[$i];
    if ($line=~m/^ATOM.{17}$chain\s{0,3}(\-?\d{1,4})/)
    {
    my $startnum=$1;
  #  print " c:$chain line:$line startnum found : $startnum\n" ;
      return $startnum;
      }
}
  }

  else
  {
   # print "no chain info**\n";
  my $cnt=@readpdb;
  for (my $j=0; $j<$cnt; $j++)
    {
    my $line=$readpdb[$j];
    if ($line=~m/^ATOM.{18}\s{0,3}(\-?\d{1,4})/)
      {
      my $startnum=$1;
    #  print " line: $j startnum found: $startnum\n" ;
      return $startnum;
      }
    }
   }
#print "no start residue found using default of 1\n";
return'1';
}

sub getgenseq
{
my %selectedhit=%{$_[0]};
my $range=$_[1]; 
my $Newquery=$_[2];
my $fastafile=$_[3]; 
my $seqdone=0;
my $seqid=$selectedhit{pdbid};
  
#print "FTP:$ftpdir\n";

  #print "SEQID:$seqid\n";
  my @selectedres=@{$selectedhit{resnum}};
  my $restot=@selectedres;
  my $res1= $selectedres[0] ;
  my $res2= $selectedres[$restot-1] ;
  my $seqnumber=@seqrecord;
  my $chain=$selectedhit{chain};
  #print "chain:: $chain\n";
  #print "pdb:$pdb\n"; 
 my $startresidue=1;
  #print "getting sequence for.. $seqid,$chain,$startresidue\n$fastafile\n"; 
my $endres=$startresidue;	
if ($Newquery)
{
#getnewseq($
}
else
{  
my $sequence=getgensequence($fastafile,$chain,$startresidue,$seqid);

 if($range)  
                   {
      my $lengthadjust=0;
      my $substrstart=$res1-$startresidue-$range;
		if ($substrstart<0)
			{
			$lengthadjust=$substrstart;
			$substrstart=0;
			}
      #print "res1:$res1 \- startres: $startresidue \- Range:$range \= substart:$substrstart\n";
      my $substrlength=(($res2-$res1)+1)+(2*$range)-$lengthadjust;
      #print "\(res2:$res2 \- res1 $res1\) \+ \(2xrange:$range\) \= $substrlength\n";
      $sequence=substr($sequence,$substrstart,$substrlength);
	 #print  "pdb start res:$startresidue-$endres, firstresidue:$res1, substring start:$substrstart\n";
                    }
     return "$sequence\n";
}
}
sub getgensequence
{
my $fastafile=$_[0];
my $chain=$_[1];
my $startresidue=$_[2];
my $seqid=$_[3];

my $foundid=0;
my $sequence="";
open(FASTAFILE,"<$fastafile");
@fastafile=<FASTAFILE>;
close FASTAFILE;

$linecnt=@fastafile;

for(my $i=0;$i<$linecnt;$i++)
{

my $line=$fastafile[$i];
if ($line=~m/>$seqid/ && $foundid==0)
	{
	$foundid=1;
	#print "$seqid found\n";
	}
elsif((!($line=~m/>.*/)) && $foundid==1)
	{
	$sequence=$sequence.$line;
	}
elsif($line=~m/>.*/ && $foundid==1)
        {
	#print "SEQ:$sequence\n";
        return $sequence;
        }
else
	{}
}
return $sequence;
}
sub getsequence
{
my $openpdb=$_[0];
my $chain=$_[1];
my $res1=$_[2];
my $res2=$_[3];
my $try=$res1;
  open(PDBFILE, $openpdb);
  my @readpdb = <PDBFILE>;
  close PDBFILE;
my $error="";
#print "getsequence:OPENING:$openpdb,\nCHAIN:$chain,RES:$res1,RES:$res2,\n";
 my $seq="";
 my $lsttry=0;
my $startflag=0;

if ($chain && $chain ne "-" && $chain ne " ")
  {
	#print "***chain true != \"-\" && != \" \"\n";
    my $cnt=@readpdb;
 OUTLOOP1: for (my $i=0;$i<$cnt;$i++)
    {
    my $line=$readpdb[$i];
	#if ($chain eq "H")
	#{
       # print "$openpdb\n$line $try\n$seq\n";
	#}
	#ATOM      9  N   GLN A   3      29.002  -0.982   7.755  1.00 61.18           N
	
	#$splittry=$try;
        #$alttry="";
       # @splittry=split(//,$splittry);
	#shift(@splittry);
	##for $each(@splittry)
	#{
	#$alttry=$alttry.$eac;
	#}

    if ($line=~m/^ATOM.{13}(\w{3})\s$chain\s{0,3}$try\s/)
	{
	$error=$error."A adding residue $1: $line";
	my $add=$1." ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	$startflag=1;
	}
    elsif ($line=~m/^ATOM.{13}(\w{3})\s{2,5}$try\s/)
	{
	$error=$error."B adding residue $1: $line\n";
	my $add=$1." ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	$startflag=1;
	}
   #  elsif ($line=~m/^ATOM.{13}(\w{3})\s$chain\s{0,4}$alttry\w\s/)
	#{
	#print "C adding residue $1: $line";
	#my $add=$1." ";
    	#$seq=$seq.$add;
	#$endres++;
	#$lsttry=$alttry;
	#$try++;
	#$startflag=1;
	#}
    
    elsif ($startflag==1&&$line=~m/^ATOM.{13}(\w{3})\s$chain\s{0,4}$lsttry/)
	{
	$error=$error."D try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	}
    elsif ($startflag==1&&$line=~m/^ATOM.{13}(\w{3})\s$chain\s{0,3}(\d{0,5})/ && $try<$2)
	{
	$error=$error."E try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	my $add="UNK ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	$i--;
	}
    elsif ($startflag==1&&$line=~m/^ATOM.{13}(\w{3})\s\w\s{0,4}/)
	{
	$error=$error."F try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	}
   elsif($line=~m/^ATOM/)
	{
	$error=$error."G error line does not match any options\n$openpdb\n$line\nTRY:$try\nSEQ:$seq\nCHAIN:$chain\n";
  	return;
	}
     elsif ($line=~m/^HETATM.{11}(\w{3})\s$chain\s{0,3}$try\s/)
	{
		#print "1try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	my $add="UNK ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	}
     elsif ($line=~m/^HETATM.{11}(\w{3})\s{2,5}$try\s/)
	{
		$error=$error."H try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	my $add="UNK ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	}
   
      elsif ($line=~m/^HETATM/)
	{
		$error=$error."I try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	}
      elsif ($startflag==1&&$line=~m/^ATOM.{15}\w\s.\s.*/)
	{
	$error=$error."J try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	}
      elsif ($startflag==1&&$line=~m/^TER/)#||$line=~m/^TER.*/
	{
	$error=$error."K try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
	if ($seq=~m/\s*\n/)
	{
	$seq="terminated on line:$line\n";
	#print "end of molecule\n";
	}
	last OUTLOOP1;
	}
     elsif ($line=~m/^TERM/||$line=~m/^ENDMDL/||$line=~m/^END/)#||$line=~m/^TER.*/
	{
	$error=$error."L try=$try lsttry=$lsttry adding residue UNK: \n$line\n";
#	print "end of molecule\n";
	last OUTLOOP1;
	}
    
    }
  }
elsif (!$chain||$chain eq "-"||$chain eq " ")
  {
	#print "*****chain = \"-\" || = \" \"\n";
   # print "no chain info\n";
  my$cnt=@readpdb;
  OUTLOOP2:for (my $i=0;$i<$cnt;$i++)
    {
    my $line=$readpdb[$i];
        if ($line=~m/^ATOM.{13}(\w{3})\s[\w\s]\s{0,3}$try\s/)
	{
	#print "$i:$1:\n$line";
	#print "ATOMXXXXXXXXXXXXXWWW W$try:$lsttry\n";
	#print "adding residue $1: $line\n";
	my $add=$1." ";
    	$seq=$seq.$add;
	$endres++;
        $lsttry=$try;
	$try++;
	$startflag=1;
	}

    elsif ($startflag==1&&$line=~m/^ATOM.{13}(\w{3})\s[\w\s]\s{0,3}$lsttry\s/)
	{
	#print "$i:-:\n$line";
	#print "ATOMXXXXXXXXXXXXXWWW W$try:$lsttry\n";
	}
    elsif ($startflag==1&&$line=~m/^ATOM.{13}(\w{3})\s[\w\s]\s{0,3}/)
	{
	#print "try=$try lsttry=$lsttry adding residue UNK: $line\n";
	my $add="UNK ";
    	$seq=$seq.$add;
	$endres++;
	#$lsttry=$try;
	$try++;
	$i--;
	}
     elsif ($line=~m/^HETATM.{11}(\w{3})\s[\w\s]\s{0,3}$try\s/)
	{
	#print "try=$try lsttry=$lsttry adding residue UNK: $line\n";
	my $add="UNK ";
    	$seq=$seq.$add;
	$endres++;
	$lsttry=$try;
	$try++;
	}
    elsif ($line=~m/^HETATM/)
	{}
     elsif ($line=~m/^TERM/||$line=~m/^ENDMDL/||$line=~m/^TER/||$line=~m/^END/)
	{last OUTLOOP2}
    }
  }
if ($seq eq "")
{
return $error;
}
 #print "$openpdb\n$seq\n";
return $seq;
}

sub sortmotifdata
{
my @hitlist=@{$_[0]};
my $count=@hitlist;
my @outlist;

#print "Removeing hits with a motif not starting with D...\n";
for (my$i=0;$i<$count;$i++)
{
if(($i%100)==0)
{
my $percent=(100/$count)*$i;
$percent=int($percent);
#print "$percent\% complete\n"
}
my %hit=%{$hitlist[$i]};
my $motif=$hit{motif};
my $r1=substr($motif,0,1);
#print "first motif residue:$r1\n";
if($r1 eq 'D')
{
push (@outlist,\%hit);
}
else
{
print"$hit{file} removed motif is: $motif\n";
print "seqence is:$hit{allseq}.\n";
}
}
return @outlist;
}

sub addmotifdata
{
my @hitlist=@{$_[0]};
my $count=@hitlist;
my @outlist;
my $ftpdir=$_[1];


#print "adding motif sequence data...\n";
for (my$i=0;$i<$count;$i++)
{
if(($i%100)==0)
{
my $percent=(100/$count)*$i;
$percent=int($percent);
print "$percent\% complete\n"
}
my %hit=%{$hitlist[$i]};
my $motif=getseq(\%hit,1,"",$ftpdir);
my $allseq=getseq(\%hit,0,"",$ftpdir);
#print "SEQ:$motif";
my $motif=substr($motif,1,5);
#print "MOTIF:$motif\n";
$hit{motif}=$motif;
$hit{allseq}=$allseq;
push (@outlist,\%hit);
}

return @outlist;

}

sub modclans
{
my @scoplist=@{$_[1]};
my @hitlist=@{$_[0]};
my $scopcnt=@scoplist;
my $hitcnt=@hitlist;
my @outarray;

for (my $i=0; $i<$scopcnt; $i++)
	{
	my %scop=%{$scoplist[$i]};
	my $scoptest=$scop{scop};
	
	for (my $j=0; $j<$hitcnt; $j++)
		{
		my %hit=%{$hitlist[$j]};
		my $scophit=$hit{scop};
		if ($scoptest eq $scophit)
			{
			push(@outarray, \%hit);
			}
		}
	}
return @outarray;

}

sub sortnegatives2
{
my @scoplist=@{$_[1]};
my @hitlist=@{$_[0]};
my $scopcnt=@scoplist;
my $hitcnt=@hitlist;
my @outarray;



for (my $i=0; $i<$scopcnt; $i++)
	{
	my %scop=%{$scoplist[$i]};
	my $scoptest=$scop{scop};
	#print "searching hitlist for scop:$scoptest\n";
	if ($i%100==1)
	{
	#print "$i/$scopcnt\n";
	}
	my $scopflag=0;
	for (my $j=0; $j<$hitcnt; $j++)
		{
		my %hit=%{$hitlist[$j]};
		my $scophit=$hit{scop};
		#print "searching against scop:$scophit\n";
		if ($scoptest eq $scophit)
			{
			#print "***$scoptest = $scophit removed from list***\n";
			$scopflag=1;
			}
		}
	if (!$scopflag)
		{
		#print "scop $scophit not found adding to list\n";
		push (@outarray, \%scop)
		}
	}
#print @outarray;
return @outarray;

}
sub sortnegatives1
{
my @scoplist=@{$_[1]};
my @hitlist=@{$_[0]};
my $scopcnt=@scoplist;
my $hitcnt=@hitlist;
my @outarray;



for (my $i=0; $i<$scopcnt; $i++)
	{
	my %scop=%{$scoplist[$i]};
	$scoptest=$scop{pdbid};
	#print "searching hitlist for scop:$scoptest\n";
	if ($i%100==1)
	{
	#print "$i/$scopcnt\n";
	}
	my $scopflag=0;
	for (my $j=0; $j<$hitcnt; $j++)
		{
		my %hit=%{$hitlist[$j]};
		$scophit=$hit{pdbid};
		#print "searching against scop:$scophit\n";
		if ($scoptest eq $scophit)
			{
			#print "***$scoptest = $scophit removed from list***\n";
			$scopflag=1;
			}
		}
	if (!$scopflag)
		{
		#print "scop $scophit not found adding to list\n";
		push (@outarray, \%scop)
		}
	}
#print @outarray;
return @outarray;

}
sub sortpositives
{
my @scoplist=@{$_[1]};#hand made scop list
my @hitlist=@{$_[0]};#positive hit list
my $scopcnt=@scoplist;
my $hitcnt=@hitlist;
my @outarray;

for (my $i=0; $i<$scopcnt; $i++)
	{
	my %scop=%{$scoplist[$i]};
	my $scoptest=$scop{scop};
	#$scoptest=~m/(\w\.\d{1,3})\.\d{1,3}/;
	#$scoptest=$1;
	print "searching hitlist for scop:$scoptest\n";
	my $scopflag=0;
	for (my $j=0; $j<$hitcnt; $j++)
		{
		my %hit=%{$hitlist[$j]};
		my $scophit=$hit{scop};
		$scophit=~m/(\w\.\d{1,3}\.\d{1,3})\.\d{1,3}/;
		#$scophit=~m/(\w\.\d{1,3})\.\d{1,3}\.\d{1,3}/;
		$scophit=$1;
		if ($scoptest eq $scophit)
			{
			print "***$scoptest = $scophit added to list***\n";
			$scopflag=1;
			push (@outarray, \%hit)
			}
		}
	if (!$scopflag)
		{
		print "scop $scophit not found removing from list\n";
		
		}
	}
#print @outarray;
return @outarray;

}
sub filtersuperfamily

{
my @input=@{$_[0]};
my $filter=$_[1];
my @output;
my $incnt=@input;
my $addcnt=0;

for (my $i=0; $i<$incnt; $i++)
	{
	my %scopI=%{$input[$i]};
	my $scopI=$scopI{scop};
	$scopI=~m/(\w\.\d{1,3}\.\d{1,3})\.\d{1,3}/;
	my $scopI2=$1;
	print "**$scopI2\n"; 
	my $outcnt=@output;

	print "searching hitlist for scop:$scopI\n";
	my $scopflag=0;
	for (my $j=0; $j<$outcnt; $j++)
		{
		my %scopO=%{$output[$j]};
		my $scopO=$scopO{scop};
		$scopO=~m/(\w\.\d{1,3}\.\d{1,3})\.\d{1,3}/;
		my $scopO2=$1;

		if ($scopI2 eq $scopO2)
			{
			print "***$scopI2 = $scopO2 removed from list***\n";
			$scopflag=1;
			}
		}
	if (!$scopflag)
		{
		if($addcnt % $filter==1)
			{
			print "Filter superfamily scop $scopI2 not found adding to list\n";
			push (@output, \%scopI);
			}
		$addcnt++;
		}
	}
return @output;

}
sub getseq
{
my %selectedhit=%{$_[0]};
my $range=$_[1]; 
my $Newquery=$_[2];
my $ftpdir=$_[3]; 
my $seqdone=0;
my $pdb=lc($selectedhit{pdbid});
  
#print "FTP:$ftpdir\n";

  print "$pdb\n";
  my @selectedres=@{$selectedhit{resnum}};
  my $restot= @selectedres;
  my $res1= $selectedres[0] ;
  my $res2= $selectedres[$restot-1] ;
  my $seqnumber=@seqrecord;
  my $chain=$selectedhit{chain};
  #print "chain:: $chain\n";
 # print "pdb:$pdb\n";
print "getseq:";
pdbftp($ftpdir,$pdb);
  my $openpdb ="$ftpdir$pdb.pdb";
 
 my $startresidue=getstartresidue($openpdb,$chain) ;
  #print "getseq:getting sequence for.. $pdb,$chain,$startresidue\n$openpdb\n"; 
my $endres=$startresidue;	
if ($Newquery)
{
#getnewseq($
}
else
{  
my $selectedchain=getsequence($openpdb,$chain,$startresidue);
#print "seq:$selectedchain\n";
 my $oneaacode=iub3to1($selectedchain);
  #print $oneaacode;
 
  if($range)  
                   {
      my $lengthadjust=0;
      my $substrstart=$res1-$startresidue-$range;
		if ($substrstart<0)
			{
			$lengthadjust=$substrstart;
			$substrstart=0;
			}
      #print "res1:$res1 \- startres: $startresidue \- Range:$range \= substart:$substrstart\n";
      my $substrlength=(($res2-$res1)+1)+(2*$range)-$lengthadjust;
      #print "\(res2:$res2 \- res1 $res1\) \+ \(2xrange:$range\) \= $substrlength\n";
      $oneaacode=substr($oneaacode,$substrstart,$substrlength);
	 #print  "pdb start res:$startresidue-$endres, firstresidue:$res1, substring start:$substrstart\n";
                    }

     return "$oneaacode\n";
}
}

1;
