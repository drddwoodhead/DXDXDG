sub saveDTdata
{
my @results1=@{$_[0]};
#my @results2=@{$_[1]};
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
my $avesize=$_[18];
my $avehyd=$_[19];
my $resultsdir=$_[20];
my $runident=$_[21];

my $filename=$resultsdir.$runident.".data";
my $top=50;
my $bottom=5;

#print ("file created: $filename\n");
open(RESULTS,">$filename"); # open new file
saveDTset(\@results1,$aasize,$aatype,$ss,$con,$sol,$abstdata,$iddata,$runident,$aa,$size,$abssdata,$hyd,$resn,\@aasizearray,\@sizearray,$aahyd,\@aahydarray,$avesize,$avehyd);
print "file created:\n$filename\n";
}

sub saveDTset
{
my @results1=@{$_[0]};
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
my $avesize=$_[18];
my $avehyd=$_[19];

my $res1cnt=@results1;
#my $reslcnt=7;

for (my $i=0; $i<$res1cnt; $i++)
	{
	
  	my %each=%{$results1[$i]};
  	my @pdblist=@{$each{pdblist}};
  	my $scop=$each{scop};
	my $restotc=@{$each{resnum}};
        my $pdbiddata="ID:".$each{'repid'};
	my $bvalue=$each{bvalue};
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
	my $avesizestr;
	my $avehydstr;
	

	my @resdata=@{$each{resdata}};
	my $resdatacnt=@resdata;
		#print "resdata: @resdata\n";
		for (my $i=0;$i<$resdatacnt;$i++)
			{
			my %res=%{$resdata[$i]};
			#my @res=%res;
			$pdbiddata="ID:".$res{'id'};
			$polar=$polar.int($res{'perP'}+0.5).",";
			#print "***$polar\n";
			$npolar=$npolar.int($res{'perN'}+0.5).",";
			$neg=$neg.int($res{'perNe'}+0.5).",";
			$pos=$pos.int($res{'perPl'}+0.5).",";
			$gly=$gly.int($res{'perG'}+0.5).",";
			$avesizestr=$avesizestr.int($res{'meansize'}+0.5).",";
			}
	if($aasize eq "v")
	{
	my $cnt=@aasizearray;
	for (my $i=0;$i<$cnt;$i++)
				{
				my $aathreshold=$aasizearray[$i];
				my @resdata=@{$each{'resdata'.$aathreshold}};
				my %res=%{$resdata[$i]};
				#my @res=%res;
				if ($aathreshold != 0)
					{
					$vsmall=$vsmall.int($res{'perS'}+0.5).",";
					$vlarge=$vlarge.int($res{'perL'}+0.5).",";
					print "vsmall:$vsmall,vlarge:$vlarge\n";
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
			#my @res=%res;
			$pdbiddata="ID:".$res{'id'};
			$small=$small.int($res{'perS'}+0.5).",";
			$large=$large.int($res{'perL'}+0.5).",";
			}
	}
	
	if($aahyd eq "v")
	{
	#print "***aahyd=V\n***";
	for (my $i=0;$i<13;$i++)
				{
				my $aathreshold=$aahydarray[$i];
				my @resdata=@{$each{'resdatah'.$aathreshold}};
				my %res=%{$resdatah[$i]};
				#my @res=%res;
				if ($aathreshold != 0)
					{
					$vpho=$vpho.int($res{'perPO'}+0.5).",";
					$vphil=$vphil.int($res{'perPI'}+0.5).",";
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
			my %res=%{$resdata[$i]};
			#my @res=%res;
			$pho=$pho.int($res{'perPO'}+0.5).",";
			$phil=$phil.int($res{'perPI'}+0.5).",";
			$avehydstr=$avehydstr.sprintf('%.1f',$res{'meanhyd'}+0.05).",";
			#print "perPO:$res{'perPO'}, perPI:$res{'perPI'}\n";
			}
	}
		
	for (my $i=0;$i<13;$i++)
				{
				my @sizedata=@{$each{sizedata}};
				my %res=%{$sizedata[$i]};
				#my @res=%res;
				if ($sizearray[$i]==1)
				{$vsizestr=$vsizestr.$res{'size'}.",";}
				}
	if ($resn==0)	
			{
			my @aasizedata=@{$each{sizedata}};
			$sizedatacnt=@sizedata;
			for (my $i=0;$i<$aasizedatacnt;$i++)
				{
				my %res=%{$aasizedata[$i]};
				#my @res=%res;
				$sizestr=$sizestr.$res{'size'}.",";
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
				$sizestr=$sizestr.$res{'size'}.",";
				#print "!Binds OUTSIZE:".$res{'size'}.",".$res{'type'}."\n :$sizestr:\n";
				
				
			}
	
	my @hyddata=@{$each{hyddata}};
		$hyddatacnt=@hyddata;
	for (my $i=0;$i<$aahyddatacnt;$i++)
			{
			my %res=%{$hyddata[$i]};
			#my @res=%res;
			$hydstr=$hydstr.$res{'hyd'}.",";
			#print "!Binds OUTSIZE:".$res{'size'}.",".$res{'type'}."\n :$sizestr:\n";
			}	


	my @AAdata=@{$each{AAdata}};
		$AAdatacnt=@AAdata;
	for (my $i=0;$i<$AAdatacnt;$i++)
			{
			my %res=%{$AAdata[$i]};
			#my @res=%res;
			$aastr=$aastr.$res{'AA'}.",";
			#print "OUTAA:".$res{'size'}.",".$res{'type'}."\n";
			}

	my @absdata=@{$each{absdata}};
		$absdatacnt=@absdata;
	for (my $i=0;$i<$absdatacnt;$i++)
			{
			my %res=%{$absdata[$i]};
			#my @res=%res;
			$abssize=$abssize.$res{'size'}.",";
			$abstype=$abstype.$res{'type'}.","; 
			#print "OUTABS:".$res{'size'}.",".$res{'type'}."\n";
			}
	my %ssdata=%{$each{SSdata}};
	#$SSD=$SSD.$ssdata{updist}.",".$ssdata{downdist}.",";
	$SSS=$SSS.$ssdata{uplength}.",".$ssdata{downlength}.",".$ssdata{ssscore}.",".$ssdata{updist}.",".$ssdata{downdist}.",".$ssdata{uptype}.",".$ssdata{downtype}.",";
	my @sable=@{$each{sable}};
	my $sablecnt=@sable;
	for (my $i=0; $i<$sablecnt; $i++)
		{
		%sableres=%{$sable[$i]};
		$solacc=$solacc.$sableres{acc}.",";
		}

	my @conserved=@{$each{conserved}};
	 $conservedtest0=testconserved(\@conserved,0,0).",";
	 $conservedtest1=testconserved(\@conserved,1,1).",";
	 $conservedtest2=testconserved(\@conserved,2,2).",";
	 $conservedtest3=testconserved(\@conserved,3,3).",";
	 $conservedtest4=testconserved(\@conserved,4,4).",";
	 $conservedtest5=testconserved(\@conserved,5,5).",";
	 $conservedtest6=testconserved(\@conserved,6,6).",";
	 $conservedtest7=testconserved(\@conserved,7,7).",";
	 $conservedtest8=testconserved(\@conserved,8,8).",";
	 $conservedtest9=testconserved(\@conserved,9,9).",";
	 $conservedtest10=testconserved(\@conserved,10,10).",";
	 $conservedtest11=testconserved(\@conserved,20,11).",";
	 $conservedtest12=testconserved(\@conserved,40,21).",";
	 $conservedtest13=testconserved(\@conserved,60,41).",";
	 $conservedtest14=testconserved(\@conserved,80,61).",";
	 $conservedtest15=testconserved(\@conserved,100,81).",";
	
	if($iddata==0){$pdbiddata="";}
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
	if($avesize==0){$avesizestr="";}
	if($avehyd==0){$avehydstr="";}

	my $printline=$pdbiddata.$abssize.$abstype.$small.$large.$vsmall.$vlarge.$polar.$npolar.$neg.$pos.$gly.$conservedtest0.$conservedtest1.$conservedtest2.$conservedtest3.$conservedtest4.$conservedtest5.$conservedtest6.$conservedtest7.$conservedtest8.$conservedtest9.$conservedtest10.$conservedtest11.$conservedtest12.$conservedtest13.$conservedtest14.$conservedtest15.$SSD.$SSS.$solacc.$sizestr.$vsizestr.$aastr.$hydstr.$pho.$phil.$vpho.$vphil.$avesizestr.$avehydstr;
	#print $printline; 
	$printline=$printline."$bvalue\n";
	print RESULTS $printline; 
	#print "::$printline";
  	}
}

sub saveDTnames
{
my $aasize=$_[0];
my $aatype=$_[1];
my $ss=$_[2];
my $con=$_[3];
my $sol=$_[4];
my $abstdata=$_[5];
my $iddata=$_[6];
my $runident=$_[7];
my $aa=$_[8];
my $size=$_[9];
my $abssdata=$_[10];
my $hyd=$_[11];
my $res=$_[12];
my @aasizearray=@{$_[13]};
my @sizearray=@{$_[14]};
my $aahyd=$_[15];
my @aahydarray=@{$_[16]};
my $avesize=$_[17];
my $avehyd=$_[18];
my $resultsdir=$_[19];
my $runident=$_[20];

my $filename=$resultsdir.$runident.".names";
my $top=50;
my $bottom=5;

open(RESULTS,">$filename"); # open new file
print RESULTS "B+,B-.\n\n";
if($iddata==1)
	{}
if($abssdata==1)
	{
	print RESULTS"ABS-4: L, S, X.\nABS-3: L, S, X.\nABS-2: L, S, X.\nABS-1: L, S, X.\nABS0: L, S, X.\nABS1: L, S, X.\nABS2: L, S, X.\nABS3: L, S, X.\nABS4: L, S, X.\nABS5: L, S, X.\nABS6: L, S, X.\nABS7: L, S, X.\nABS8: L, S, X.\n";
	}
if($abstdata==1)
	{
	print RESULTS"ABT-4: P,N,+,-,G, X.\nABT-3: P,N,+,-,G, X.\nABT-2: P,N,+,-,G, X.\nABT-1: P,N,+,-,G, X.\nABT0: P,N,+,-,G, X.\nABT1: P,N,+,-,G, X.\nABT2: P,N,+,-,G, X.\nABT3: P,N,+,-,G, X.\nABT4: P,N,+,-,G, X.\nABT5: P,N,+,-,G, X.\nABT6: P,N,+,-,G, X.\nABT7: P,N,+,-,G, X.\nABT8: P,N,+,-,G, X.\n";
	}
if($aasize!=0)
	{	
	print RESULTS "AST-4S: continuous.\nAST-3S: continuous.\nAST-2S: continuous.\nAST-1S: continuous.\nAST0S: continuous.\nAST1S: continuous.\nAST2S: continuous.\nAST3S: continuous.\nAST4S: continuous.\nAST5S: continuous.\nAST6S: continuous.\nAST7S: continuous.\nAST8S: continuous.\nAST-4L: continuous.\nAST-3L: continuous.\nAST-2L: continuous.\nAST-1L: continuous.\nAST0L: continuous.\nAST1L: continuous.\nAST2L: continuous.\nAST3L: continuous.\nAST4L: continuous.\nAST5L: continuous.\nAST6L: continuous.\nAST7L: continuous.\nAST8L: continuous.\n";
	}
if ($aasize eq "v")
{
my $loopv=@aasizearray;
for (my $i=0;$i<$loopv;$i++)
	{
	$each=$aasizearray[$i];
	if($each!=0)
		{
		my $res=$i-4;
		print RESULTS "AST$res"."S".": continuous.\n";
		}
	}
for (my $i=0;$i<$loopv;$i++)
	{
	$each=$aasizearray[$i];
	if($each!=0)
		{
		my $res=$i-4;
		print RESULTS "AST$res"."L".": continuous.\n";
		}	
	}
}


if($aatype==1)
	{
	print RESULTS "ATG-4P: continuous.\nATG-3P: continuous.\nATG-2P: continuous.\nATG-1P: continuous.\nATG0P: continuous.\nATG1P: continuous.\nATG2P: continuous.\nATG3P: continuous.\nATG4P: continuous.\nATG5P: continuous.\nATG6P: continuous.\nATG7P: continuous.\nATG8P: continuous.\nATG-4N: continuous.\nATG-3N: continuous.\nATG-2N: continuous.\nATG-1N: continuous.\nATG0N: continuous.\nATG1N: continuous.\nATG2N: continuous.\nATG3N: continuous.\nATG4N: continuous.\nATG5N: continuous.\nATG6N: continuous.\nATG7N: continuous.\nATG8N: continuous.\nATG-4Ne: continuous.\nATG-3Ne: continuous.\nATG-2Ne: continuous.\nATG-1Ne: continuous.\nATG0Ne: continuous.\nATG1Ne: continuous.\nATG2Ne: continuous.\nATG3Ne: continuous.\nATG4Ne: continuous.\nATG5Ne: continuous.\nATG6Ne: continuous.\nATG7Ne: continuous.\nATG8Ne: continuous.\nATG-4Pl: continuous.\nATG-3Pl: continuous.\nATG-2Pl: continuous.\nATG-1Pl: continuous.\nATG0Pl: continuous.\nATG1Pl: continuous.\nATG2Pl: continuous.\nATG3Pl: continuous.\nATG4Pl: continuous.\nATG5Pl: continuous.\nATG6Pl: continuous.\nATG7Pl: continuous.\nATG8Pl: continuous.\nATG-4G: continuous.\nATG-3G: continuous.\nATG-2G: continuous.\nATG-1G: continuous.\nATG0G: continuous.\nATG1G: continuous.\nATG2G: continuous.\nATG3G: continuous.\nATG4G: continuous.\nATG5G: continuous.\nATG6G: continuous.\nATG7G: continuous.\nATG8G: continuous.\n";
	}
if($con==1)
	{
	print RESULTS"CON0: True,False.\nCON1: True,False.\nCON2: True,False.\nCON3: True,False.\nCON4: True,False.\nCON5: True,False.\nCON6: True,False.\nCON7: True,False.\nCON8: True,False.\nCON9: True,False.\nCON10: True,False.\nCON11: True,False.\nCON12: True,False.\nCON13: True,False.\nCON14: True,False.\nCON15: True,False.\n";}

if($ss==1)
	{
	print RESULTS"SSDu: continuous.\nSSDd: continuous.\nSSSm: continuous.\nSSLu: continuous.\nSSLd: continuous.\nSSTu:H,E,-.\nSSTd:H,E,-.\n";
	}
if($sol==1)
	{
	print RESULTS"SOL-4: continuous.\nSOL-3: continuous.\nSOL-2: continuous.\nSOL-1: continuous.\nSOL0: continuous.\nSOL1: continuous.\nSOL2: continuous.\nSOL3: continuous.\nSOL4: continuous.\nSOL5: continuous.\nSOL6: continuous.\nSOL7: continuous.\nSOL8: continuous.\n";
	}
if($size==1 && $res==0)
	{
	print RESULTS"ABW-4: continuous.\nABW-3: continuous.\nABW-2: continuous.\nABW-1: continuous.\nABW0: continuous.\nABW1: continuous.\nABW2: continuous.\nABW3: continuous.\nABW4: continuous.\nABW5: continuous.\nABW6: continuous.\nABW7: continuous.\nABW8: continuous.\n";
	}
if($size==1 && $res!=0)
	{	
	my $aa=$res-5;
	print RESULTS "ABW$aa: continuous.\n";
	}
if ($size eq "v")
	{
	my $loopv=@sizearray;
for (my $i=0;$i<$loopv;$i++)
	{
	$each=$sizearray[$i];
	if($each!=0)
		{
		my $res=$i-4;
		print RESULTS "ABW$res:$each: continuous.\n";
		}
	}
	}
if($avesize==1)
	{
	print RESULTS"AveABW-4: continuous.\nAveABW-3: continuous.\nAveABW-2: continuous.\nAveABW-1: continuous.\nAveABW0: continuous.\nAveABW1: continuous.\nAveABW2: continuous.\nAveABW3: continuous.\nAveABW4: continuous.\nAveABW5: continuous.\nAveABW6: continuous.\nAveABW7: continuous.\nAveABW8: continuous.\n";
	}	

if($aa==1)
	{
	print RESULTS"ABR-4: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR-3: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR-2: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR-1: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR0: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR1: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR2: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR3: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR4: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR5: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR6: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR7: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\nABR8: G,A,S,P,V,T,C,L,I,N,D,Q,K,E,M,H,F,R,Y,W,X.\n";
	}
if($hyd==1)
	{
	print RESULTS"ABH -4: continuous.\nABH -3: continuous.\nABH -2: continuous.\nABH -1: continuous.\nABH 0: continuous.\nABH 1: continuous.\nABH 2: continuous.\nABH 3: continuous.\nABH 4: continuous.\nABH 5: continuous.\nABH 6: continuous.\nABH 7: continuous.\nABH 8: continuous.\n";
	}
if($avehyd==1)
	{
	print RESULTS"HYA -4: continuous.\nHYA -3: continuous.\nHYA -2: continuous.\nHYA -1: continuous.\nHYA 0: continuous.\nHYA 1: continuous.\nHYA 2: continuous.\nHYA 3: continuous.\nHYA 4: continuous.\nHYA 5: continuous.\nHYA 6: continuous.\nHYA 7: continuous.\nHYA 8: continuous.\n";
	}
if($aahyd!=0)
	{	
	print RESULTS "HYCpho-4: continuous.\nHYCpho-3: continuous.\nHYCpho-2: continuous.\nHYCpho-1: continuous.\nHYCpho0: continuous.\nHYCpho1: continuous.\nHYCpho2: continuous.\nHYCpho3: continuous.\nHYCpho4: continuous.\nHYCpho5: continuous.\nHYCpho6: continuous.\nHYCpho7: continuous.\nHYCpho8: continuous.\n\HYCphi-4: continuous.\n\HYCphi-3: continuous.\n\HYCphi-2: continuous.\n\HYCphi-1: continuous.\n\HYCphi0: continuous.\n\HYCphi1: continuous.\n\HYCphi2: continuous.\n\HYCphi3: continuous.\n\HYCphi4: continuous.\n\HYCphi5: continuous.\n\HYCphi6: continuous.\n\HYCphi7: continuous.\n\HYCphi8: continuous.\n";
	}
if ($aahyd eq "v")
{
my $loopv=@aahydarray;
for (my $i=0;$i<$loopv;$i++)
	{
	$each=$aahydarray[$i];
	if($each!=0)
		{
		my $res=$i-4;
		print RESULTS "HYCpho$res: continuous.\n";
		}
	}
for (my $i=0;$i<$loopv;$i++)
	{
	$each=$aahydarray[$i];
	if($each!=0)
		{
		my $res=$i-4;
		print RESULTS "HYCphi$res: continuous.\n";
		}	
	}
}
print "file created:\n$filename\n";
}

sub writetreedata
{
#my $dtfile=$_[1];
my %errordata=%{$_[0]};
my $treefile=$_[1];
my $skip=$_[2];
open (TREEFILE,">>$treefile");
print "$treefile\n";
#open (DTFILE,">>$dtfile");
my @tree=@{$errordata{tree}};
#my $per=$errordata{pererror};
#my $ernum=$errordata{numerror};
my $dtcnt=@tree;
#print TREEFILE "$per,$ernum,";
#print DTFILE "ERRORS\n";
print TREEFILE "$skip\t";
for ($i=0;$i<$dtcnt;$i++)
{
my $line=$tree[$i];
print TREEFILE "$line\t";
}
print TREEFILE "\n";

}

sub resid
{
my $number=$_[0]+2;

if ($number<10)
{
return "$number";
}
else
{
return $number;
}
}
sub getDTerrors
{
my $filename=$_[0];
my @familydata=@{$_[1]};
my %dtdata=%{$_[2]};
my %dtdata2=%{$_[3]};
my $runident=$_[4];
my $dtresultsdir=$_[5];
my $skip=$_[6];


my $errorfile="$dtresultsdir$runident.errors";
my $treefile="$dtresultsdir$runident.trees";
my $statsfile="$dtresultsdir$runident.stats";
my $mergefile="$dtresultsdir$runident.merge";
#print "OPEN: $filename\n";
open(FILE,"<$filename");
my @file = <FILE>;
close FILE;
my @errors;
my @tree;
my $linecnt=@file;
my $treeflag=0;
#print "Lines=$linecnt\n";
for (my $i=0;$i<$linecnt; $i++)
{
$line=$file[$i];
if ($treeflag==1)
	{
	if($line=~m/^Tree saved.*/||$line=~m/^Simplified Decision Tree:.*/)
		{
		$treeflag=0;
		}
	else
		{
		chomp($line);
		push(@tree,$line);
		}
}
elsif($line=~m/^<error>(.*)<\\error>/)
{
print "Error: $line ::$1\n";
my $error=$1;
$error=~s/\s//g;
push(@errors, $error);
}
elsif($line=~m/^\s*\d*\s*(\d*)\((\s*\d*\.\d)%\).*<</)
{
#	  34    3( 5.2%)     13   10(17.2%)    (32.9%)   <<
print "$line\n";
$ernum=$1;
$per=$2;
}
elsif($line=~m/^Decision Tree:.*/)
{
$treeflag=1;
}

else 
{
#print "no match:: $line";
}
}

#print "**skip:$skip\n";
my @errorout=geterrordata(\@errors,\@familydata);
my %errordata;
%errordata=(errors,\@errorout,pererror,$per,numerror,$ernum,tree,\@tree);
writeerrordata(\%errordata,$filename,$errorfile,$skip);
writetreedata(\%errordata, $treefile,$skip);
writestatsdata(\%dtdata,$statsfile,$skip);
writestatsdata(\%dtdata2,$statsfile,$skip);


return %errordata;
}



sub geterrordata
{
my @errors=@{$_[0]};
my @familydata=@{$_[1]};
my $errorcnt=@errors;
my @errorout;

for (my $i=0; $i<$errorcnt;$i++)
{
my $errornum=$errors[$i];
my %hit=%{$familydata[$errornum]};
$hit{'dtid'}=$errornum;
push (@errorout, \%hit); 
}
return @errorout;

}

1;

