sub saveSVMdata
{
my @results=@{$_[0]};
my $aasize=$_[1];
my $aatype=$_[2];
my $ss=$_[3];
my $con=$_[4];
my $sol=$_[5];
my $abstdata=$_[6];
my $iddata=$_[7];
my $runident=$_[8];
my $aa=$_[9];
my $size=$_[10];
my $abssdata=$_[11];
my $hyd=$_[12];
my $resn=$_[13];
my @aasizearray=@{$_[14]};
my @sizearray=@{$_[15]};
my $aahyd=$_[16];
my @aahydarray=@{$_[17]};
my @aatypearray=@{$_[18]};
my $resdir=$_[19];
my $runid=$_[20];


my $top=50;
my $bottom=5;

#print ("file created: $filename\n");
saveSVMset(\@results,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runident,$aa,$size,$abssdata,$hyd,$resn,\@aasizearray,\@sizearray,$aahyd, \@aahydarray,\@aatypearray,$resdir,$runid);


}

sub saveSVMset
{
my @results=@{$_[0]};
my $aasize=$_[1];
my $aatype=$_[2];
my $ss=$_[3];
my $con=$_[4];
my $sol=$_[5];
my $abstdata=$_[6];
my $iddata=$_[7];
my $runident=$_[8];
my $aa=$_[9];
my $size=$_[10];
my $abssdata=$_[11];
my $hyd=$_[12];
my $resn=$_[13];
my @aasizearray=@{$_[14]};
my @sizearray=@{$_[15]};
my $aahyd=$_[16];
my @aahydarray=@{$_[17]};
my @aatypearray=@{$_[18]};
my $resdir=$_[19];
my $runid=$_[20];

my $res1cnt=@results;
#my $reslcnt=7;
my $filename=$resdir.$runid.".data";
open(RESULTS,">$filename"); # open new file
print "#file created:\n$filename\n";
for (my $i=0; $i<$res1cnt; $i++)
	{
	
  	my %each=%{$results[$i]};
  	my @pdblist=@{$each{pdblist}};
  	my $scop=$each{scop};
	my $restotc=@{$each{resnum}};
        my $pdbiddata="\#ID:".$each{'repid'};
	my $bindvalue=$each{bvalue};
	my $set;
	if ($bindvalue eq "B+")
		{
		$set="1 ";
		print "SETB+:$bindvalue:$set\n";
		}
	elsif ($bindvalue eq "B-")
		{
		$set="-1 ";
		print "SETB-:$bindvalue:$set\n";
		}
	my $abssize;
	my $abstype;
  	my $small;
	#$small="S";
	my $large;
	#$large="L";
	my $polar;
	#$polar="P";
	my $npolar;
	#$npolar="N";
	my $neg;
	#$neg="-";
	my $pos;
	#$pos="+";
	my $gly;
	#$gly="G";
	my $Conserved="False,";
	my $SSD;
	#$SSD="D";
	my $SSS;
	#$SSS="^";
	my $solacc;
	#$solacc="A";
	my $sizestr;
	my $aastr;
	my $hydstr;
	my $vsmall;
	my $vlarge;
	my $vsizestr;
	my $pho;
	my $phil;
	my $vpho;
	my $vphil;

	my @resdata=@{$each{resdata}};
	
	if ($aatype eq "v")
	{
	for (my $i=0;$i<17;$i++)
				{
				my $aathreshold=$aatypearray[$i];
				my @resdata=@{$each{'resdata'}};
				my %res=%{$resdata[$i]};
				#my @res=%res;
				if ($aathreshold != 0)
					{
					my $resid=resid($i);
					my $id="01";
					#$pdbiddata="\#ID:XXX".$res{'id'}."\n";
					$polar=$polar.$id."00".$resid.":".int($res{'perP'}+0.5)." ";
					$npolar=$npolar.$id."01".$resid.":".int($res{'perN'}+0.5)." ";
					$neg=$neg.$id."02".$resid.":".int($res{'perNe'}+0.5)." ";
					$pos=$pos.$id."03".$resid.":".int($res{'perPl'}+0.5)." ";
					$gly=$gly.$id."04".$resid.":".int($res{'perG'}+0.5)." "; 
					}
					
				}
	}
	
	else
	{
	my $resdatacnt=@resdata;
		#print "resdata: @resdata\n";
		#print "resdatacnt: $resdatacnt\n";
		for (my $i=0;$i<$resdatacnt;$i++)
			{
			my %res=%{$resdata[$i]};
			#my @res=%res;
			#print "@res/n";
			my $resid=resid($i);
			my $id="01";
			#$pdbiddata="\#ID:XXX".$res{'id'}."\n";
			$polar=$polar.$id."00".$resid.":".int($res{'perP'}+0.5)." ";
			$npolar=$npolar.$id."01".$resid.":".int($res{'perN'}+0.5)." ";
			$neg=$neg.$id."02".$resid.":".int($res{'perNe'}+0.5)." ";
			$pos=$pos.$id."03".$resid.":".int($res{'perPl'}+0.5)." ";
			$gly=$gly.$id."04".$resid.":".int($res{'perG'}+0.5)." ";
			}
	}
	if($aasize eq "v")
	{
	
	for (my $i=0;$i<17;$i++)
				{
				my $aathreshold=$aasizearray[$i];
				my @resdata=@{$each{'resdata'.$aathreshold}};
				my %res=%{$resdata[$i]};
				#my @res=%res;
				if ($aathreshold != 0)
					{
					my $resid=resid($i);
					my $id="00";
					
					$vsmall=$vsmall.$id."00".$resid.":".int($res{'perS'}+0.5)." ";
					$vlarge=$vlarge.$id."01".$resid.":".int($res{'perL'}+0.5)." ";
					my $test= $res{'perS'};
					#print "::TESTAASV::$test ::\n"; 
					}
					
				}	
	}
	else
	{
	my @resdata=@{$each{resdata}};
	my $resdatacnt=@resdata;
		#print "resdata: @resdata\n";
		for (my $i=0;$i<$resdatacnt;$i++)
			{
			my %res=%{$resdata[$i]};
			my $resid=resid($i);
			my $id="00";
			my @res=%res;
			#print "RES:@res\n";
			#$pdbiddata="ID:".$res{'id'};
			my $test=$res{'perL'};
			my $test2=$res{'perS'};
			$small=$small.$id."00".$resid.":".int($test2+0.5)." ";
			$large=$large.$id."01".$resid.":".int($test+0.5)." ";
			
			#print "::TESTAAS$aasize ::L$test ::S$test2\n"; 
			#print "$small\n$large\n";
			}
	}
	
	if($aahyd eq "v")
	{
	#print "***aahyd=V\n***";
	for (my $i=0;$i<17;$i++)
				{
				my $aathreshold=$aahydarray[$i];
				my @resdata=@{$each{'resdatah'.$aathreshold}};
				my %res=%{$resdatah[$i]};
				#my @res=%res;
				if ($aathreshold != 0)
					{
					my $resid=resid($i);
					my $id="12";
					
					
					$vpho=$vpho.$id."00".$resid.":".int($res{'perPO'}+0.5)." ";
					$vphil=$vphil.$id."01".$resid.":".int($res{'perPI'}+0.5)." ";
					}
					
				}	
	}
	else
	{
	my @resdata=@{$each{resdatah}};
	my $resdatacnt=@resdata;
		#print "resdata: @resdata\n";
		for (my $i=0;$i<$resdatacnt;$i++)
			{
			my $resid=resid($i);
			my $id="12";
			my %res=%{$resdata[$i]};
			#my @res=%res;
			$pho=$pho.$id."00".$resid.":".int($res{'perPO'}+0.5)." ";
			$phil=$phil.$id."01".$resid.":".int($res{'perPI'}+0.5)." ";
			}
	}
		
	for (my $i=0;$i<17;$i++)
				{
				my @sizedata=@{$each{sizedata}};
				my %res=%{$sizedata[$i]};
				#my @res=%res;
				if ($sizearray[$i]==1)
				{
				my $resid=resid($i);
					my $id="03";
				$vsizestr=$vsizestr.$id."00".$resid.":".$res{'size'}." ";}
				}
	if ($resn==0)	
			{
			my @aasizedata=@{$each{sizedata}};
			$sizedatacnt=@sizedata;
			for (my $i=0;$i<$aasizedatacnt;$i++)
				{
				my $resid=resid($i);
				my $id="03";
				my %res=%{$aasizedata[$i]};
				#my @res=%res;
				$sizestr=$sizestr.$id."00".$resid.":".$res{'size'}." ";
				#print "!Binds OUTSIZE:".$res{'size'}.",".$res{'type'}."\n :$sizestr:\n";
				}
			}
	elsif	($resn!=0)
			{
			my @aasizedata=@{$each{sizedata}};
			$aasizedatacnt=@aasizedata;	
				$myres=$resn-1;
				my %res=%{$aasizedata[$myres]};
				#my @res=%res;
				my $resid=resid($i);
				my $id="03";
				$sizestr=$sizestr.$id."00".$resid.":".$res{'size'}." ";
				#print "!Binds OUTSIZE:".$res{'size'}.",".$res{'type'}."\n :$sizestr:\n";
				
				
			}
	
	my @hyddata=@{$each{hyddata}};
		$hyddatacnt=@hyddata;
	for (my $i=0;$i<$aahyddatacnt;$i++) #ABH DATA ABsolute Hydrophobicity
			{
			my $resid=resid($i);
			my $id="07";
			my %res=%{$hyddata[$i]};
			#my @res=%res;
			$hydstr=$hydstr.$id."00".$resid.":".$res{'hyd'}." ";
			#print "!Binds OUTSIZE:".$res{'size'}.",".$res{'type'}."\n :$sizestr:\n";
			}	


	my @AAdata=@{$each{AAdata}};
		$AAdatacnt=@AAdata;
	for (my $i=0;$i<$AAdatacnt;$i++)
			{
			my %res=%{$AAdata[$i]};
			#my @res=%res;
			my $resid=resid($i);
			my $id="02".$resid;
			my $aaid=aaid($res{'AA'});
			$aastr=$aastr.$resid.$aaid.":1 ";
			#print "OUTAA:".$res{'size'}.",".$res{'type'}."\n";
			}

	my @absdata=@{$each{absdata}};
		$absdatacnt=@absdata;
	for (my $i=0;$i<$absdatacnt;$i++)
			{
			my %res=%{$absdata[$i]};
			#my @res=%res;
			my $resid=resid($i);
			my $id="04".$resid;
			if ($res{'size'} eq 'small')
				{	
				$sid='00';
				$abssize=$abssize.$sid.$id.":1 ";
				}
			else
				{}
			
			}
	for (my $i=0;$i<$absdatacnt;$i++)
			{
			my %res=%{$absdata[$i]};
			#my @res=%res;
			my $resid=resid($i);
			my $id="04".$resid;
			if ($res{'size'} eq 'large')
				{	
				$sid='01';
				$abssize=$abssize.$sid.$id.":1 ";
				}
			else
				{}
			}
	
	for (my $i=0;$i<$absdatacnt;$i++)
			{
			my $resid=resid($i);
			$id="05".$resid;
			my $tid=typeid($res{'type'});
			$abstype=$abstype.$tid.$id.":1 "; 
			#print "OUTABS:".$res{'size'}.",".$res{'type'}."\n";
			}
	my %ssdata=%{$each{SSdata}};
	#$SSD=$SSD.$ssdata{updist}.",".$ssdata{downdist}.",";
	$SSS=$SSS."110102:".$ssdata{uplength}." "."110202:".$ssdata{downlength}." "."110300:".$ssdata{ssscore}." "."110101:".$ssdata{updist}." "."110201:".$ssdata{downdist}." "."110103:".gettype($ssdata{uptype})." "."110203:".gettype($ssdata{downtype})." ";
	
	my @sable=@{$each{sable}};
	my $sablecnt=@sable;
	#print "SABLECNT:$sablecnt\n";
	for (my $i=0; $i<$sablecnt; $i++)
		{
		%sableres=%{$sable[$i]};
		my $resid=resid($i);
		$id="05";
		my $test = $sableres{acc};
		#print "SABLERES:$test\n";
		$solacc=$solacc.$id."00".$resid.":".$sableres{acc}." ";
		}

	my @conserved=@{$each{conserved}};
	 $conservedtest0="100000:".testconserved2(\@conserved,0,0)." ";
	 $conservedtest1="100100:".testconserved2(\@conserved,1,1)." ";
	 $conservedtest2="100200:".testconserved2(\@conserved,2,2)." ";
	 $conservedtest3="100300:".testconserved2(\@conserved,3,3)." ";
	 $conservedtest4="100400:".testconserved2(\@conserved,4,4)." ";
	 $conservedtest5="100500:".testconserved2(\@conserved,5,5)." ";
	 $conservedtest6="100600:".testconserved2(\@conserved,6,6)." ";
	 $conservedtest7="100700:".testconserved2(\@conserved,7,7)." ";
	 $conservedtest8="100800:".testconserved2(\@conserved,8,8)." ";
	 $conservedtest9="100900:".testconserved2(\@conserved,9,9)." ";
	 $conservedtest10="101000:".testconserved2(\@conserved,10,10)." ";
	 $conservedtest11="101100:".testconserved2(\@conserved,20,11)." ";
	 $conservedtest12="101200:".testconserved2(\@conserved,40,21)." ";
	 $conservedtest13="101300:".testconserved2(\@conserved,60,41)." ";
	 $conservedtest14="101400:".testconserved2(\@conserved,80,61)." ";
	 $conservedtest15="101500:".testconserved2(\@conserved,100,81)." ";
	
	#if($iddata==0){$pdbiddata="";}
	if($abstdata==0){$abstype="";}
	if($abssdata==0){$abssize="";}
	if($aasize==0||$aasize eq "v"){$small="";$large="";}
	if($aasize ne "v"){$vsmall="";$vlarge="";}
	if($aatype==0){$polar="";$npolar="";$neg="";$pos="";$gly="";}
	if($ss==0){$SSD="";$SSS="";}
	if($con==0){$conservedtest0="",$conservedtest1="";$conservedtest2="";$conservedtest3="";$conservedtest4="";$conservedtest5="";$conservedtest6="";$conservedtest7="";$conservedtest8="";$conservedtest9="";$conservedtest10="";$conservedtest11="";$conservedtest12="";$conservedtest13="";$conservedtest14="";$conservedtest15="";}
	if($sol==0){$solacc="";}
	if($aa==0){$aastr="";}
	if($size==0){$sizestr="";}
	if($size ne "v"){$vsizestr="";}
	if($hyd==0){$hydstr="";}
	if($aahyd==0||$aahyd eq "v"){$pho="";$phil="";}
	if($aahyd ne "v"){$vpho="";$vphil="";}
	my $printline=$abssize.$abstype.$small.$large.$vsmall.$vlarge.$polar.$npolar.$neg.$pos.$gly.$conservedtest0.$conservedtest1.$conservedtest2.$conservedtest3.$conservedtest4.$conservedtest5.$conservedtest6.$conservedtest7.$conservedtest8.$conservedtest9.$conservedtest10.$conservedtest11.$conservedtest12.$conservedtest13.$conservedtest14.$conservedtest15.$SSD.$SSS.$solacc.$sizestr.$vsizestr.$aastr.$hydstr.$pho.$phil.$vpho.$vphil;
	#print $printline;
	$printline=sortdata($printline); 
	#print "**SORT$printline\n";
	if (!$set)
	{
	print RESULTS "0 ";
	}
	$printline=$set.$printline." ".$pdbiddata;
	print RESULTS $printline."\n"; 
	print "::$printline\n";
  	}
#@testresults=<RESULTS>;
#print "RESULTSFILE:@testresults\n";
}

sub resid
{
my$id=$_[0];
$id=$id+2;
if ($id <10 )
{
$id="0".$id;
return $id;
}
else
{
return $id;
}
}

sub aaid
{
my $aa=$_[0];
$residue=uc($residue);

if($residue eq 'G') 
{
return "00";
}
elsif($residue eq 'A') 
{
return "01";
}
elsif($residue eq 'S') 
{
return "02";
}
elsif($residue eq 'P') 
{
return "03";
}
elsif($residue eq 'V')
{
return "04";
}
elsif($residue eq 'T')
{
return "05";
}
elsif($residue eq 'C')
{
return "06";
}
elsif($residue eq 'L')
{
return "07";
}
elsif($residue eq 'I')
{
return "08";
}
elsif($residue eq 'N')
{
return "09";
}
elsif($residue eq 'D')
{
return "10";
}
elsif($residue eq 'Q')
{
return "11";
}
elsif($residue eq 'K')
{
return "12";
}
elsif($residue eq 'E')
{
return "13";
}
elsif($residue eq 'M')
{
return "14";
}
elsif($residue eq 'H')
{
return "15";
}
elsif($residue eq 'F')
{
return "16";
}
elsif($residue eq 'R')
{
return "17";
}
elsif($residue eq 'Y')
{
return "18";
}
elsif($residue eq 'W')
{
return "19";
}
elsif($residue eq "X"|| $residue eq ""||$residue eq "-"||$residue eq " ") 
{
return "";
}
else
{
#print "**$residue**\n";
}
}

sub typeid
{
my $aa=$_[0];
$residue=uc($residue);

if($residue eq '+') 
{
return "00";
}
elsif($residue eq '-') 
{
return "01";
}
elsif($residue eq 'P') 
{
return "02";
}
elsif($residue eq 'N') 
{
return "03";
}
elsif($residue eq 'G')
{
return "04";
}

}

sub sortdata
{
my $data=$_[0];

my @data=split(/ /, $data);
@data = sort { $a <=> $b } @data;

my $list= "@data";
return $list;

}

sub writevectordata
{
my $vectorfile=$_[1];
#my $dtfile=$_[1];
my %errordata=%{$_[0]};
open (VECTORFILE,">>$vectorfile");

#open (DTFILE,">>$dtfile");
my @vector=@{$errordata{vectors}};
#my $per=$errordata{pererror};
#my $ernum=$errordata{numerror};
my $vcnt=@vector;
#print TREEFILE "$per,$ernum,";
#print DTFILE "ERRORS\n";
for ($i=0;$i<$vcnt;$i++)
{
my $line=$vector[$i];
print VECTORFILE "$line\t";
}
print VECTORFILE "\n";

}

sub getSVMerrors
{
my $filename2=$_[0];
my $filename1=$_[1];
my @familydata=@{$_[2]};
my %svmdata=%{$_[3]};
my %svmdata2=%{$_[4]};
my $runident=$_[5];
my $svmresultsdir=$_[6];
my $skip=$_[7];


my $errorfile="$svmresultsdir$runident.errors";
my $vectorfile="$svmresultsdir$runident.vec";
my $statsfile="$svmresultsdir$runident.stats";
my $mergefile="$svmresultsdir$runident.merge";





open(FILE1,"<$filename1");
my @file1 = <FILE1>;
close FILE1;
my @vector;
my $linecnt1=@file1;
my $vflag=0;
#print "Lines=$linecnt1\n";
for (my $i=0;$i<$linecnt1; $i++)
{
$line=$file1[$i];
if ($vflag==1)
	{
	if($line=~m/^/)
		{
		$vflag=0;
		}
	else
		{
		chomp($line);
		push(@vector,$line);
		}
}
elsif($line=~m/^\d*\.\d* # threshold b, each following line is a SV \(starting with alpha\*y\)\*/)
{
$vflag=1;
}


else 
{
#print "no match:: $line";
}
}

open(FILE2,"<$filename2");
my @file2 = <FILE2>;
close FILE2;
my $linecnt2=@file2;
my %info;
my $miss;
my $errorest;
my $recallest;
my $precisionest;
my @errors;

#print "Lines=$linecnt\n";
for (my $i=0;$i<$linecnt2; $i++)
{
$line=$file2[$i];
#print "$line";
if($line=~m/^Optimization finished \((\d*) misclassified, maxdiff=\d*\.\d*\)./)
{
#print "miss :: $1\n";
$miss=$1;
}
elsif($line=~m/^XiAlpha-estimate of the error: error<=(\d*\.\d*)\% \(rho=\d*.\d*,depth=\d*\)/)
{
#	  34    3( 5.2%)     13   10(17.2%)    (32.9%)   <<
#print "XIerror :: $1\n";
$errorest=$1;
}
elsif($line=~m/^XiAlpha-estimate of the recall: recall=>(\d*\.\d*)\% \(rho=\d*\.\d*,depth=\d*\)/)
{
#	  34    3( 5.2%)     13   10(17.2%)    (32.9%)   <<
#print "XI recall :: $1\n";
$recallest=$1;
}
elsif($line=~m/^XiAlpha-estimate of the precision: precision=>(\d*\.\d*)\% \(rho=\d*\.\d*,depth=\d*\)/)
{
#	  34    3( 5.2%)     13   10(17.2%)    (32.9%)   <<
#print "XI precision :: $1\n";
$precisionest=$1;
}
elsif($line=~m/^<error>(.*)<\\error>/)
{
#print "Error: $line ::$1\n";
my $error=$1;
$error=~s/\s//g;
$error=~s/ID://;
push(@errors, $error);
}
else 
{
#print "no match:: $line";
}
}



#my @errorout=geterrordata(\@errors,\@familydata);
#my %errordata;
%svmerrordata=(errors,\@errors,miss,$miss,errorest,$errorest,recallest,$recallest,precisionest,$precisionest);
writesvmerrordata(\%svmerrordata,$filename,$errorfile);
writestatsdata(\%svmdata,$statsfile,$skip);
writestatsdata(\%svmdata2,$statsfile,$skip);
mergefile($errorfile,$vectorfile,$statsfile,$mergefile);

return %svmerordata;
}

sub writesvmerrordata
{
my $file=$_[1];
my %errordata=%{$_[0]};
my $errorfile=$_[2];
open (ERRORFILE,">>$errorfile");

open (FILE,">>$file");
my @errors=@{$errordata{errors}};
my $miss=$errordata{miss};
my $errorest=$errordata{errorest};
my $recallest=$errordata{recallest};
my $precisionest=$errordata{precisionest};
my $ercnt=@errors;
print ERRORFILE "$miss\t$errorest\t$recallest\t$precisionest\t";
for (my $i=0;$i<$ercnt;$i++)
{
$err=$errors[$i];
print ERRORFILE "$err\t";
}
#print DTFILE "ERRORS\n";
#for ($i=0;$i<$ercnt;$i++)
{
my %hit=%{$errors[$i]};
my $dtid=$hit{dtid};
my $id=$hit{repid};
my $scop=$hit{scop};
print FILE "$dtid\t$id\t$scop\t";
print ERRORFILE "$dtid\t$id\t$scop\t";

}
print ERRORFILE "\n";

}

sub getsupportvectors
{
my $file= $_[0];
my $attribfile=$_[1];
my @vectorarray;

my $modelfile="$file.model";
#print "openfile:$modelfile\n";

open (MODELFILE, "$modelfile");

my @modelfile=<MODELFILE>;
my $linecnt=@modelfile;

#print "@modelfile\n";

for (my$i=0;$i<$linecnt;$i++)
{
my %supportvector;
my $line=$modelfile[$i];
#print "*$line";
if ($line=~m/^(\d{1,4})\s\#\snumber\sof\ssupport\svectors\splus\s1/)
{
$vectorcnt=$1;
#print "VECTORCOUNT= $vectorcnt\n";
}
else
{

}
#print "LINE:$line";
if ($line=~m/^-{0,1}\d{1,4}\.\d{1,35}e{0,1}-{0,1}\d{0,4}\s.*/)
{
$line=~m/^(-{0,1}\d{1,4}\.\d{1,35}e{0,1}-{0,1}\d{0,4})\s.*\s#(.*)/;
$alphay=$1;
$id=$2;
$line=~s/^.*?\s//; 
$line=~s/#.*//;
%supportvector=(alphay,$alphay,id,$id);
#print "alpha*y=$alphay, ID=$id\n";
@vectors=split(/\s/,$line);
#print "**VECTORS @vectors\n";
$veccnt=@vectors;

for (my $i=0; $i<$veccnt; $i++)
	{
	my $vecpnt=$vectors[$i];
	$vecpnt=~m/(\d{1,6}):(\d{1,6})/;
	my $attribid=$1;
	my $attribvalue=$2;
	$supportvector{$attribid}=$attribvalue;
	}
my @supportvector=%supportvector;
#print "**EXVEC HASH:@supportvector\n";
my %vec=%supportvector;
@vec=%vec;

push (@vectorarray,\%vec);
}
else
{#print "NO VECTORS FOUND\n";
}
}
my %vectors= %{$vectorarray[1]};
my @vectors=%vectors;
#print "exampleHASHfromarray: @vectors\n";
writesupportvectorfile(\@vectorarray,$file,$attribfile);

#return @vectorarray;
}

sub writesupportvectorfile
{
my @supportvectors=@{$_[0]};
my $file=$_[1];
my $vectab=$_[2];

my %vectors= %{$supportvectors[1]};
my @vectors=%vectors;
#print "exampleHASH2fromarray: @vectors\n";


my @attribs=getsvattribs(\@supportvectors);

my %attriblist=getsvattriblist(\@supportvectors,\@attribs);

#if(-e $vectab)
#{
#print "AFF\n";
#}
#else
#{
#print "NO ATTRIB FILE FOUND\n";
open (MODELTAB, ">>$vectab");
my $attribcnt=@attribs;
print "attribcnt:$attribcnt\n";
for (my $i=0;$i<$attribcnt;$i++)
{
my $attribid=getattribid($attribs[$i]);
#print "ATT:$attribs[$i]:,ATTRIBID:$attribid:\n";
print MODELTAB "$attribid\t";
}
print MODELTAB "\n";
#}
#open (MODELTAB,">>$vectab");

my $attribcnt=@attribs;

for (my $i=0;$attribcnt>$i;$i++)
	{
	my $attrib=$attribs[$i];
	my $result=$attriblist{$attrib};
	#print "AAT:$attrib,RE:$result\n";
	print MODELTAB "$result\t";
	}
print MODELTAB "\n";
}

sub getsvattriblist
{
my @supportvectors=@{$_[0]};
my @attribs=@{$_[1]};

my %vectors= %{$supportvectors[1]};
my @vectors=%vectors;
#print "exampleHASH4fromarray: @vectors\n";


my $attribcnt=@attribs;

my %attribhash;

for (my $i=0;$i<$attribcnt;$i++)
{
my $att=$attribs[$i];
#print "ATT:$att\n";

my $veccnt=@supportvectors;
#print "veccnt:$veccnt\n";
$attribhash{$att}=0;
for (my $j=0;$j<$veccnt;$j++)
		{
		my %vector=%{$supportvectors[$j]};
		if ($vector{$att}>10)
			{$attribhash{$att}++;}
		}


} 
@attribhash=%attribhash;
#print "ATTHASH @attribhash\n";
return %attribhash;
}


sub getsvattribs
{
my @vectorarray=@{$_[0]};
my %vectors= %{$vectorarray[1]};
my @vectors=%vectors;
my @outlist;


my %vectors= %{$vectorarray[1]};
my @vectors=%vectors;
#print "exampleHASH3fromarray: @vectors\n";

delete $vectors{id};
delete $vectors{alphay};
my @vectors=%vectors;
#print "VECTORS:@vectors\n";
my $veccnt=@vectors;

for (my $i=0; $veccnt>$i; $i=$i+2)
{
push(@outlist, $vectors[$i]);
}
#print "VECTORATTRIBS:@outlist\n";
return @outlist;
}

sub getattribid
{

my $query=$_[0];
if ($query==2)
{
return "AST-4S";
}
elsif ($query==3)
{
return "AST-3S";
}
elsif ($query==4)
{
return "AST-2S";
}
elsif ($query==5)
{
return "AST-1S";
}
elsif ($query==6)
{
return "AST0S";
}
elsif ($query==7)
{
return "AST1S";
}
elsif ($query==8)
{
return "AST2S";
}
elsif ($query==9)
{
return "AST3S";
}
elsif ($query==10)
{
return "AST4S";
}
elsif ($query==11)
{
return "AST5S";
}
elsif ($query==12)
{
return "AST6S";
}
elsif ($query==13)
{
return "AST7S";
}
elsif ($query==14)
{
return "AST8S";
}

elsif ($query==102)
{
return "AST-4L";
}
elsif ($query==103)
{
return "AST-3L";
}
elsif ($query==104)
{
return "AST-2L";
}
elsif ($query==105)
{
return "AST-1L";
}
elsif ($query==106)
{
return "AST-0L";
}
elsif ($query==107)
{
return "AST1L";
}
elsif ($query==108)
{
return "AST2L";
}
elsif ($query==109)
{
return "AST3L";
}
elsif ($query==110)
{
return "AST4L";
}
elsif ($query==111)
{
return "AST5L";
}
elsif ($query==112)
{
return "AST6L";
}
elsif ($query==113)
{
return "AST7L";
}
elsif ($query==114)
{
return "AST8L";
}

elsif ($query==10002)
{
return "ATG-4P";
}
elsif ($query==10003)
{
return "ATG-3P";
}
elsif ($query==10004)
{
return "ATG-2P";
}
elsif ($query==10005)
{
return "ATG-1P";
}
elsif ($query==10006)
{
return "ATG0P";
}
elsif ($query==10007)
{
return "ATG1P";
}
elsif ($query==10008)
{
return "ATG2P";
}
elsif ($query==10009)
{
return "ATG3P";
}
elsif ($query==10010)
{
return "ATG4P";
}
elsif ($query==10011)
{
return "ATG5P";
}
elsif ($query==10012)
{
return "ATG6P";
}
elsif ($query==10013)
{
return "ATG7P";
}
elsif ($query==10014)
{
return "ATG8P";
}


elsif ($query==10102)
{
return "ATG-4N";
}
elsif ($query==10103)
{
return "ATG-3N";
}
elsif ($query==10104)
{
return "ATG-2N";
}
elsif ($query==10105)
{
return "ATG-1N";
}
elsif ($query==10106)
{
return "ATG-0N";
}
elsif ($query==10107)
{
return "ATG1N";
}
elsif ($query==10108)
{
return "ATG2N";
}
elsif ($query==10109)
{
return "ATG3N";
}
elsif ($query==10110)
{
return "ATG4N";
}
elsif ($query==10111)
{
return "ATG5N";
}
elsif ($query==10112)
{
return "ATG6N";
}
elsif ($query==10113)
{
return "ATG7N";
}
elsif ($query==10114)
{
return "ATG8N";
}

elsif ($query==10202)
{
return "ATG-4Ne";
}
elsif ($query==10203)
{
return "ATG-3Ne";
}
elsif ($query==10204)
{
return "ATG-2Ne";
}
elsif ($query==10205)
{
return "ATG-1Ne";
}
elsif ($query==10206)
{
return "ATG-0Ne";
}
elsif ($query==10207)
{
return "ATG1Ne";
}
elsif ($query==10208)
{
return "ATG2Ne";
}
elsif ($query==10209)
{
return "ATG3Ne";
}
elsif ($query==10210)
{
return "ATG4Ne";
}
elsif ($query==10211)
{
return "ATG5Ne";
}
elsif ($query==10212)
{
return "ATG6Ne";
}
elsif ($query==10213)
{
return "ATG7Ne";
}
elsif ($query==10214)
{
return "ATG8Ne";
}


elsif ($query==10302)
{
return "ATG-4Pl";
}
elsif ($query==10303)
{
return "ATG-3Pl";
}
elsif ($query==10304)
{
return "ATG-2Pl";
}
elsif ($query==10305)
{
return "ATG-1Pl";
}
elsif ($query==10306)
{
return "ATG-0Pl";
}
elsif ($query==10307)
{
return "ATG1Pl";
}
elsif ($query==10308)
{
return "ATG2Pl";
}
elsif ($query==10309)
{
return "ATG3Pl";
}
elsif ($query==10310)
{
return "ATG4Pl";
}
elsif ($query==10311)
{
return "ATG5Pl";
}
elsif ($query==10312)
{
return "ATG6Pl";
}
elsif ($query==10313)
{
return "ATG7Pl";
}
elsif ($query==10314)
{
return "ATG8Pl";
}

elsif ($query==10402)
{
return "ATG-4G";
}
elsif ($query==10403)
{
return "ATG-3G";
}
elsif ($query==10404)
{
return "ATG-2G";
}
elsif ($query==10405)
{
return "ATG-1G";
}
elsif ($query==10406)
{
return "ATG-0G";
}
elsif ($query==10407)
{
return "ATG1G";
}
elsif ($query==10408)
{
return "ATG2G";
}
elsif ($query==10409)
{
return "ATG3G";
}
elsif ($query==10410)
{
return "ATG4G";
}
elsif ($query==10411)
{
return "ATG5G";
}
elsif ($query==10412)
{
return "ATG6G";
}
elsif ($query==10413)
{
return "ATG7G";
}
elsif ($query==10414)
{
return "ATG8G";
}

elsif ($query==30002)
{
return "AAW-4";
}
elsif ($query==30003)
{
return "AAW-3";
}
elsif ($query==30004)
{
return "AAW-2";
}
elsif ($query==30005)
{
return "AAW-1";
}
elsif ($query==30006)
{
return "AAW0";
}
elsif ($query==30007)
{
return "AAW1";
}
elsif ($query==30008)
{
return "AAW2";
}
elsif ($query==30009)
{
return "AAW3";
}
elsif ($query==30010)
{
return "AAW4";
}
elsif ($query==30011)
{
return "AAW5";
}
elsif ($query==30012)
{
return "AAW6";
}
elsif ($query==30013)
{
return "AAW7";
}
elsif ($query==30014)
{
return "AAW8";
}


elsif ($query==110102)
{
return "SSLu";
}
elsif ($query==110202)
{
return "SSLd";
}
elsif ($query==110101)
{
return "SSDu";
}
elsif ($query==110201)
{
return "SSDd";
}
elsif ($query==110103)
{
return "SSTu";
}
elsif ($query==110203)
{
return "SSTd";
}
elsif ($query==110300)
{
return "SSSm";
}


elsif ($query==100100)
{
return "CON1";
}
elsif ($query==100200)
{
return "CON2";
}
elsif ($query==100300)
{
return "CON3";
}
elsif ($query==100400)
{
return "CON4";
}
elsif ($query==100500)
{
return "CON5";
}
elsif ($query==100600)
{
return "CON6";
}
elsif ($query==100700)
{
return "CON7";
}
elsif ($query==100800)
{
return "CON8";
}
elsif ($query==100900)
{
return "CON9";
}
elsif ($query==101000)
{
return "CON10";
}
elsif ($query==101100)
{
return "CON11";
}
elsif ($query==101200)
{
return "CON12";
}
elsif ($query==101300)
{
return "CON13";
}
elsif ($query==101400)
{
return "CON14";
}
elsif ($query==101500)
{
return "CON15";
}

}
sub gettype
{
$in=$_[0];
if($in eq "H")
{return 1}
elsif ($in eq "E")
{return -1}
else
{return 0}
}
1;