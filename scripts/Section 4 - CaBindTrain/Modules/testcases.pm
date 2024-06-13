
sub readrules
{
$rulefile=$_[0];

open (RULEFILE, "$rulefile");

my @rulefile=<RULEFILE>;

#print "reading rule file:$rulefile\n";
#print @rulefile;

my $filecnt=@rulefile;
my $inrule=0;
my @rules;

for (my $i=0;$i<$filecnt;$i++)
	{
	my $line=$rulefile[$i];
	if ($inrule==1)
	{
	
	my $result;
		if ($line=~m/^\s*->\s*class (.*)\s*\[.*\]/)
			{
			#print"MATCH:$line";
			$result=$1;
			$inrule=0;
			my %rule=('rule',$rule,'result',$result);
			#print "RULE:$rule RESULT:$result\n\n";
			push (@rules,\%rule);
			}
		else
			{
			#print "ELSE:$line";
			$line=~/^\s*(.*)\n/;
			$relement=$1;
			#print "RELEMENT:$relement\n";
			#chomp($relement);
			if($rule eq "")
			{
			#print"rule=\"\":-";
			$rule=$relement;
			print"$rule\n";
			}
			else	
			{$rule=$rule."&&".$relement;}
			print "*$rule\n";
			}
	}
	
	else
	{
	if ($line=~m/^\s*Rule\s\d.*/)
		{
		#print "INRULE:$line";
		$inrule=1;
		$rule="";
		}
	elsif ($line=~m/^\s*Default class:\s*(\w.)/)
		{
		my $result=$1;
		my %rule=('rule',"else",'result',$result);
		push (@rules,\%rule);
		}
	else {}
	}
	}
return @rules;
}

sub classifybydt
{
my $dtfileroot=$_[0];
my @familydata=@{$_[1]};
my $aasize=$_[2];
my @aasizearray=@{$_[3]};

my @truepositive;
my @truenegative;
my @falsepositive;
my @falsenegative;
my @positive;
my @negative;
my @errors;

my @rules=readrules("$dtfileroot.r");

$dtfileroot=~m/^.*\/(.+)/;
my $dtid=$1;
#print "*RULES*\n@rules\n";

my $res1cnt=@familydata;
#my $reslcnt=7;

for (my $i=0; $i<$res1cnt; $i++)
	{
	my %each=%{$familydata[$i]};
	my $id=$each{'repid'};
	#print "CASE:$id\n";
	}

my $familycnt=@familydata;

for (my $i=0;$i<$familycnt;$i++)
{
my %case=%{$familydata[$i]};
my $id=$case{'repid'};
my $kwnvalue=$case{bvalue};
#	print "CASE:$id\n";
my $testedcase=testdtcase(\@rules,\%case,$aasize,\@aasizearray);
my $caseid=$case{repid};
#print "Result: $testedcase\n";
$testedcase=~s/\s/""/g;

if ($testedcase=~m/^B\+/)
	{
#	print "$caseid:DTtest:B+\n";
	push (@positive,\%case);

	if ($kwnvalue eq "B+")
		{
#		print "TRUEPOSIVITVE\n\n";
		push (@truepositive,\%case);
		}
	elsif ($kwnvalue eq "B-")
	   	{
#		print "FALSEPOSIVITVE\n\n";
		push (@falsepositive,\%case);
		}
	else { #print "UNKNOWN\n\n";
	     }
	}

elsif ($testedcase=~m/^B\-/)
	{
#	print "$caseid:DTtest:B-\n";
	push (@negative,\%case);

	if ($kwnvalue eq "B+")
		{
#		print "FALSENEGATIVE\n\n";
		push (@falsenegative,\%case);
		}
	elsif ($kwnvalue eq "B-")
		{
#		print "TRUENEGATIVE\n\n";
		push (@truenegative,\%case);
		}
	else 	{
#		print "UNKNOWN\n\n";
		}
	}
	else
		{
#		print "$caseid:DTtest:ERROR\n\n";
		push (@errors,\%case);
		}
}
my $p=@positive;
my $n=@negative;
my $tnc=@truenegative;
my $tpc=@truepositive;
my $fnc=@falsenegative;
my $fpc=@falsepositive;

#print "*DT P:$p,N:$n,tnc:$tnc,tpc:$tpc,fnc:$fnc,fpc:$fpc\n";
my %classifications=('dtid',$dtid,'positive',\@positive,'negative',\@negative,'truepositive',\@truepositive,'truenegative',\@truenegative,'falsepositive',\@falsepositive,'falsenegative',\@falsenegative,'errors',\@errors);
return %classifications;
}

sub classifybysvm
{
my $svmmodelfile=$_[0];
my @familydata=@{$_[1]};
my $aasize=$_[2];
my $aatype=$_[3];
my $ss=$_[4];
my $con=$_[5];
my $sol=$_[6];
my $abstdata=$_[7];
my $iddata=$_[8];
my $runident=$_[9];
my $aa=$_[10];
my $size=$_[11];
my $abssdata=$_[12];
my $hyd=$_[13];
my $resn=$_[14];
my @aasizearray=@{$_[15]};
my @sizearray=@{$_[16]};
my $aahyd=$_[17];
my @aahydarray=@{$_[18]};
my @aatypearray=@{$_[19]};
my $resdir=$_[20];
my $runid=$_[21];


my @truepositive;
my @truenegative;
my @falsepositive;
my @falsenegative;
my @positive;
my @negative;
my @errors;

$svmmodelfile=~m/^.*\/(.+)\..*/;
my $svmid=$1;

my $res1cnt=@familydata;
#my $reslcnt=7;

for (my $i=0; $i<$res1cnt; $i++)
	{
	my %each=%{$familydata[$i]};
	my $id=$each{'repid'};
	#print "CASE:$id\n";
	}

my $familycnt=@familydata;

for (my $i=0;$i<$familycnt;$i++)
{
my $testedclass;
my %case=%{$familydata[$i]};
my $id=$case{'repid'};
	print "CASE:$id\n";
my $testedcase=testsvmcase($svmmodelfile,\%case,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runident,$aa,$size,$abssdata,$hyd,$res,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,\@aatypearray,$resdir,$runid);
my $caseid=$case{repid};
my $kwnvalue=$case{bvalue};

if ($line=~m/^-.+/)
	{$testedclass="B-";}
elsif ($line=~m/^[^-]/)  
	{$testedclass="B+";}
else 
	{$testedclass="ER";}

$testedclass=~s/\s/""/g;
print "Result: $testedcase : $testedclass\n";

if ($testedclass=~m/^B\+/)
	{
	print "$caseid:SVMtest:B+\n";
	push (@positive,\%case);
	print "kwnvalue is $kwnvalue\n";
	if ($kwnvalue eq "B+")
		{
		print "TRUEPOSIVITVE\n\n";
		push (@truepositive,\%case);
		}
	elsif ($kwnvalue eq "B-")
		{
		print "FALSEPOSIVITVE\n\n";
		push (@falsepositive,\%case);
		}
	else { print "UNKNOWN\n\n";}
	}

elsif ($testedclass=~m/^B\-/)
	{
	print "$caseid:SVMtest:B-\n";
	push (@negative,\%case);

	print "kwnvalue is $kwnvalue\n";
	if ($kwnvalue eq "B+")
		{
		print "FALSENEGATIVE\n\n";
		push (@falsenegative,\%case);
		}
	elsif ($kwnvalue eq "B-")
		{
		print "TRUENEGATIVE\n\n";
		push (@truenegative,\%case);
		}
	else {print "UNKNOWN\n\n";}
	}
else
	{
	print "$caseid:SVMtest:ERROR\n";
	push (@errors,\%case);
	}
}
my $p=@positive;
my $n=@negative;
my $tnc=@truenegative;
my $tpc=@truepositive;
my $fnc=@falsenegative;
my $fpc=@falsepositive;

#print "*SVM P:$p,N:$n,tnc:$tnc,tpc:$tpc,fnc:$fnc,fpc:$fpc\n";

my %classifications=('svmid',$svmid,'positive',\@positive,'negative',\@negative,'truepositive',\@truepositive,'truenegative',\@truenegative,'falsepositive',\@falsepositive,'falsenegative',\@falsenegative,'errors', \@errors);
return %classifications;
}

sub testsvmcase
{
my $modelfile=$_[0];
my %case=%{$_[1]};
my $aasize=$_[2];
my $aatype=$_[3];
my $ss=$_[4];
my $con=$_[5];
my $sol=$_[6];
my $abstdata=$_[7];
my $iddata=$_[8];
my $runident=$_[9];
my $aa=$_[10];
my $size=$_[11];
my $abssdata=$_[12];
my $hyd=$_[13];
my $resn=$_[14];
my @aasizearray=@{$_[15]};
my @sizearray=@{$_[16]};
my $aahyd=$_[17];
my @aahydarray=@{$_[18]};
my @aatypearray=@{$_[19]};
my $resdir=$_[20];
my $runid=$_[21];

my $id2=$case{'repid'};
	print "CASE*:$id2\n";

my @temparray=(\%case);

$modelfile=~m/^(.*\/).*/;
my $filename=$resdir."tempfile.data";
my $tempresult=$resdir."tempfile.result";
open(RESULTS,">$filename"); 
print "createing file.. $filename\n";
saveSVMset(\@temparray,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runident,$aa,$size,$abssdata,$hyd,$resn,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,\@aatypearray,$resdir,"tempfile");

`/home/dwoodhead/svm_light/svm_classify "$filename" "$modelfile" "$tempresult"`;

open (RESULT,"<$tempresult");
my @result=<RESULT>;
print "resultfile:@result\n";
$line=$result[0];
$line=~s/\n/""/g;
return $line;
}

sub testdtcase
{
my @rules=@{$_[0]};
my %case=%{$_[1]};
my $aasize=$_[2];
my @aasizearray=@{$_[3]};

my $id2=$case{'repid'};
	#print "CASE*:$id2\n";

my %results;

my @truepos;
my @falsepos;
my @trueneg;
my @falseneg;
my @posative;
my @negative;

my $rulecnt=@rules;

my @resdata=@{$case{resdata}};
#print "Resdata @resdata\n";

 for (my $i=0;$i<$rulecnt;$i++)
	{
	my %rule=%{$rules[$i]};
	my $result=$rule{result};
	my $thisrule=$rule{rule};
	if ($thisrule eq "else")
	{
#	print "else: $result\n";
	return $result;}

	my @ruleelements=split(/&&/, $thisrule);
	my $elecnt=@ruleelements;
	my $elecnt2=0;
#	print "Rule:$i/$rulecnt:$thisrule,Result:$result\n";
	for (my $i=0;$i<$elecnt;$i++)
			{
			my $element=$ruleelements[$i];
			$element=~m/^(.*)\s+([<,>,=]{1,2})\s+(\S+)/;
			#print "att:$1,QAL:$2,Val:$3\n";
			my $attval=getattval($1,\%case,$aasize,\@aasizearray);
			$attval=int($attval+0.5);
			#print "attribute value:$attval\n";
			my $qualifier=$2;
			my $value=$3;
			my $test=elementtruth($attval,$qualifier,$value);
#			print "$attval,$qualifier,$value:$test\n";
			if ($test)
				{
#				print "element $element :true\n";
				$elecnt2++;
				}
			else 
				{
#				print "element $element :false\n";
				}
			
			
			}
	if ($elecnt==$elecnt2)	
	{
#	print "all elements true rule true:$result  \n\n";
	return $result;
	}
	else
	{
#	print "not all elements true move to next rule\n\n";
	}
	}


}

sub elementtruth
 {
my $attvalue=$_[0];
my $qualifier=$_[1];
my $value=$_[2];

if ($value=~m/^\d+\W?\d*/)
{
	#print "*number comparison\n";
	#$value= sprintf("%.1f", $value);
	#$attvalue = sprintf("%.1f", $attvalue);
	 #print "$attvalue $qualifier $value\n";
	if($qualifier eq "=")
	{
		if($attvalue == $value)
			{return "1";}
		else
			{return "0";}
	}

	elsif($qualifier eq ">")
	{
		if($attvalue > $value)
			{
			#print "$attvalue > $value :True\n";
			return "1";	
			}
		else
			{
			#print "$attvalue > $value :False\n";
			return "0";
			}
	}
	elsif($qualifier eq "<")
	{
		if($attvalue < $value)
			{return "1";}
		else
			{return "0";}
	}
	elsif($qualifier eq ">=")
	{
		if($attvalue >= $value)
			{return "1";}
		else
			{return "0";}
	}
	elsif($qualifier eq "<=")
	{
		if($attvalue <= $value)
			{return "1";}
		else
			{return "0";}
	}
}
elsif ($value=~m/\w*/)
{
#print "*text comparison\n";
if($attvalue eq $value)
			{
			return "1";
			}
		else
			{return "0";}
}
}

sub getattval
{
my $attid=$_[0];
my %case=%{$_[1]};
my $aasize=$_[2];
my @aasizearray=@{$_[3]};

my $id3=$case{'repid'};
	print "CASE**:$id3\n";	
	print "attid:$attid\n";
if ($aasize=~m/\d/ && $attid=~m/AST(-?\d)(\w)/ )
{
print "AST\n";
my @resdata=@{$case{resdata}};
my %res=%{$resdata[$1+4]};
if ($2 eq "L")
	{return $res{perL};}
elsif ($2 eq "S")
	{return $res{perS};}
}
elsif ($aasize=~m/v/ && $attid=~m/AST(-?\d)(\w)/ )
{
print "variable AST\n";
my $res=$1+4;
my $thid=$aasizearray[$res];
print " $1 RES:$res threshold:$thid\n";
my @resdata=@{$case{'resdata'.$thid}};
my %res=%{$resdata[$res]};
if ($2 eq "L")
	{return $res{perL};}
elsif ($2 eq "S")
	{return $res{perS};}

}
elsif ($attid=~m/ATG(-?\d)(\w{1,2})/)
{
#print "ATG,$1,$2\n";
my @resdata=@{$case{resdata}};
my $resid=$1+4;
my %res=%{$resdata[$resid]};
if ($2 eq "P")
	{return $res{perP};}
elsif ($2 eq "N")
	{
	my $val=$res{'perN'};
#	print "RETURN $val\n";
	return $res{perN};
	}
elsif ($2 eq "Ne")
	{return $res{perNe};}
elsif ($2 eq "Pl")
	{return $res{perPl};}
elsif ($2 eq "G")
	{return $res{perG};}
}
elsif ($attid=~m/HYC(\w{3})(-?\d)/)
{
#print "1:$1 2:$2\n";
my @resdatah=@{$case{resdatah}};
my $resid=$2+4;
my %res=%{$resdatah[$resid]};
if ($1 eq "phi")
	{ #print "phi:$res{perPI}\n";
	return $res{'perPI'};}
elsif ($1 eq "pho")
	{#print "pho:$res{perPO}\n";
	return $res{'perPO'};}
}
elsif ($attid=~m/SOL(-?\d)/)
{
my @sable=@{$case{sable}};
my %sableres=%{$sable[$1+4]};
return $sableres{acc};
}
elsif ($attid=~m/CON(\d{1,2)/)
{

}
elsif ($attid=~m/SS([U,D,M])([S,D,L,T])/)
{
my %ssdata=%{$case{SSdata}};
if($1 eq 'U')
{
	if($2 eq 'D')
	{return $ssdata{updist}}	
	elsif($2 eq 'L')
	{return $ssdata{uplength}}
	elsif($2 eq 'T')
	{return $ssdata{uptype}}
}
elsif($1 eq 'D')
{

	if($2 eq 'D')
	{return $ssdata{downdist}}	
	elsif($2 eq 'L')
	{return $ssdata{downlength}}
	elsif($2 eq 'T')
	{return $ssdata{downtype}}
}
elsif ($1 eq 'M')
{
	return $ssdata{ssscore};
}
}
#elsif ($attid=~m/ABS.*/)
#{}
#elsif ($attid=~m/ABW.*/)
#{}
#elsif ($attid=~m/ABG.*/)
#{}
#elsif ($attid=~m/ABR.*/)
#{}
#elsif ($attid=~m/ABH.*/)
#{}
#elsif ($attid=~m/AWA.*/)
#{}
#elsif ($attid=~m/ATR.*/)
#{}
#elsif ($attid=~m/HYA.*/)
#{}
}

sub saveclassresults
{
my %results=%{$_[0]};				
my $familyfile=$_[1];
my $dataid=$_[2];
my $resultsdir=$_[3];
my $exp=$_[4];

$familyfile=~m/^.*\/(.+)\./;
$filename=$resultsdir.$1.$exp.".results";


my $id=$dtresults{dtid};
#print "creating file $filename\n";
open (OUTFILE, ">$filename");
print OUTFILE "Decision Tree Results:$id\n";

my $id=$svmresults{svmid};
#print "creating file $filename\n";
open (OUTFILE, ">$filename");
print OUTFILE "Decision Tree Results:$id\n";
printresults(\%results);
}

sub printresults
{
my %results=%{$_[0]};

my @positive=@{$results{positive}};
my @negative=@{$results{negative}};
my @truepositive=@{$results{truepositive}};
my @truenegative=@{$results{truenegative}};
my @falsepositive=@{$results{falsepositive}};
my @falsenegative=@{$results{falsenegative}};
my @errors=@{$results{errors}};

my $p=@positive;
my $n=@negative;
my $tnc=@truenegative;
my $tpc=@truepositive;
my $fnc=@falsenegative;
my $fpc=@falsepositive;
my $erc=@errors;

#print "*** P:$p,N:$n,E:$erc,tnc:$tnc,tpc:$tpc,fnc:$fnc,fpc:$fpc\n";

if ($tpc||$tnc||$fpc||$fnc)
{

print OUTFILE "True Positive	$tpc";
for (my $i=0;$i<$tpc;$i++)
	{
	%each=%{$truepositive[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
print OUTFILE "True Negative	$tnc";
for (my $i=0;$i<$tnc;$i++)
	{
	%each=%{$truenegative[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}

print OUTFILE "\n";
print OUTFILE "False Positive	$fpc";
for (my $i=0;$i<$fpc;$i++)
	{
	%each=%{$falsepositive[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
print OUTFILE "False Negative	$fnc";
for (my $i=0;$i<$fnc;$i++)
	{
	%each=%{$falsenegative[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
print OUTFILE "Errors	$erc";
for (my $i=0;$i<$erc;$i++)
	{
	%each=%{$errors[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
}
else
{
print OUTFILE "Positive	$p";
for (my $i=0;$i<$p;$i++)
	{
	%each=%{$positive[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
print OUTFILE "Negative	$n";
for (my $i=0;$i<$n;$i++)
	{
	%each=%{$negative[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
print OUTFILE "Errors	$erc";
for (my $i=0;$i<$erc;$i++)
	{
	%each=%{$errors[$i]};
	my $id = $each{repid};
	print OUTFILE "	$id";
	}
print OUTFILE "\n";
}



}

sub savegenresults
{
my %results=%{$_[0]};				
my $familyfile=$_[1];
my $dataid=$_[2];
my $resultsdir=$_[3];
my $exp=$_[4];

$familyfile=~m/^.*\/(.+)\./;
$filename=$resultsdir.$1.$exp.".gen";


my $id=$dtresults{dtid};
#print "creating file $filename\n";
open (OUTFILE, ">$filename");
print OUTFILE "Decision Tree Results:$id\n";

my $id=$svmresults{svmid};
#print "creating file $filename\n";
open (OUTFILE, ">$filename");
print OUTFILE "Decision Tree Results:$id\n";
printgenresults(\%results);
}

sub printgenresults
{
my %results=%{$_[0]};

my @positive=@{$results{positive}};
my @negative=@{$results{negative}};
my @errors=@{$results{errors}};

my $p=@positive;
my $n=@negative;
my $erc=@errors;

#print "*** P:$p,N:$n,E:$erc\n";

print OUTFILE "Positive	$p";

for (my $i=0;$i<$p;$i++)
	{
	%each=%{$positive[$i]};
	my $id = $each{repid};
	my $data = $each{gendata};
	print OUTFILE "\n$id	$data";
	my @rpsdata=@{$each{rpsdata}};
	my $rpscnt=@rpsdata;
	for (my $i=0;$i<$rpscnt;$i++)
		{
		my %rps=%{$rpsdata[$i]};
		#my $motifid=$rps{description};
		my $motifid=$rps{motifid};	
		print OUTFILE "	$motifid";
		}
	
	}
print OUTFILE "\n";
print OUTFILE "Negative	$n";
for (my $i=0;$i<$n;$i++)
	{
	%each=%{$negative[$i]};
	my $id = $each{repid};
	my $data = $each{gendata};
	print OUTFILE "\n$id	$data";
	my @rpsdata=@{$each{rpsdata}};
	my $rpscnt=@rpsdata;
	for (my $i=0;$i<$rpscnt;$i++)
		{
		my %rps=%{$rpsdata[$i]};
		#my $motifid=$rps{description};
		my $motifid=$rps{motifid};
		print OUTFILE "	$motifid";
		}
	}
print OUTFILE "\n";
print OUTFILE "Errors	$erc";
for (my $i=0;$i<$erc;$i++)
	{
	%each=%{$errors[$i]};
	my $id = $each{repid};
	my $data = $each{gendata};
	print OUTFILE "\n$id	$data";
	my @rpsdata=@{$each{rpsdata}};
	my $rpscnt=@rpsdata;
	for (my $i=0;$i<$rpscnt;$i++)
		{
		my %rps=%{$rpsdata[$i]};
		my $motifid=$rps{description};
		print OUTFILE "	$motifid";
		}
	}
print OUTFILE "\n";
}

sub writestatsdata
{
my %results=%{$_[0]};
my $statsfile=$_[1];
my $skip=$_[2];

my @positive=@{$results{positive}};
my @negative=@{$results{negative}};
my @truepositive=@{$results{truepositive}};
my @truenegative=@{$results{truenegative}};
my @falsepositive=@{$results{falsepositive}};
my @falsenegative=@{$results{falsenegative}};
my @errors=@{$results{errors}};

my $p=@positive;
my $n=@negative;
my $tnc=@truenegative;
my $tpc=@truepositive;
my $fnc=@falsenegative;
my $fpc=@falsepositive;
my $erc=@errors;

#print "*^* P:$p,N:$n,E:$erc,tnc:$tnc,tpc:$tpc,fnc:$fnc,fpc:$fpc\n";
#print "$statsfile\n";

print STATSFILE "$skip";

open(STATSFILE, ">>$statsfile");

if ($tpc||$tnc||$fpc||$fnc)
{
print STATSFILE "\t$tpc";
print STATSFILE "\t$tnc";
print STATSFILE "\t$fpc";
print STATSFILE "\t$fnc";
print STATSFILE "\t$erc";
print STATSFILE "\n";
}
else
{
print STATSFILE "\t$p";
print STATSFILE "\t$n";
print STATSFILE "\t$erc";
print STATSFILE "\n";
}

}
1;