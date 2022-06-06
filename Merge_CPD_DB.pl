#use warnings;
use Sort::Naturally 'nsort';

#die;

#input in this order
#pathbank, pubchem, hmdb, kegg, biocyc, chebi

$inPB = "PATHBANK_CPD_DB.txt";
$inPC = "PUBCHEM_CPD_DB.txt";
$inHM = "HMDB_CPD_DB.txt";
$inKG = "KEGG_CPD_DB.txt";
$inBC = "BIOCYC_CPD_DB.txt";
$inCH = "CHEBI_CPD_DB.txt";

open(INPB, $inPB)||die;
open(INPC, $inPC)||die;
open(INHM, $inHM)||die;
open(INKG, $inKG)||die;
open(INBC, $inBC)||die;
open(INCH, $inCH)||die;


print "INPUT PATHBANK\n";
$on=0;
while(<INPB>){
	if($_!~/\w/){next;}
	if($_ =~ /CPD\tFORMULA/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];
	for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];}
#	if($on%100000==0){last;}
#$on++;
}

$time=localtime;
print "INPUT HMDB on $on time $time\n";
while(<INHM>){	
	if($_!~/\w/){next;}
	if($_ =~ /^CPD\t/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];	
	for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];} 
#	if($on%100000==0){last;}
#$on++;
}


$time=localtime;
print "INPUT PUBCHEM on $on time $time\n";
while(<INPC>){
	if($_!~/\w/){next;}
	if($_ =~ /^CPD\t/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];
	for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];}
	#can use below to fill in missing info
        #if(exists($HoA{$cpd})){ 
	#	if($HoA{$cpd}[1] eq $stuff[1] && $stuff[1]=~/\w/){ 
	#		for my $i (0..12){
	#			if($HoA{$cpd}[$i]!~/\w/ && $stuff[$i]=~/\w/){ $HoA{$cpd}[$i]=$stuff[$i];}
	#}	}	}
	#else{ for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i]; }}
#	if($on%100000==0){last;}
#$on++;
}



$time=localtime;
print "INPUT KEGG on $on time $time\n";
while(<INKG>){	
	if($_!~/\w/){next;}
	if($_ =~ /^CPD\t/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];
	for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];}
#	if($on%100000==0){last;}
#$on++;
}

$time=localtime;
print "INPUT BIOCYC on $on time $time\n";
while(<INBC>){	
	if($_!~/\w/){next;}
	if($_ =~ /^CPD\t/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];
	for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];}
#	if($on%100000==0){last;}
#$on++;
}

$time=localtime;
print "INPUT CHEBI on $on time $time\n";
while(<INCH>){	
	if($_!~/\w/){next;}
	if($_ =~ /^CPD\t/i){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];
        for my $i (0..12){ $HoA{$cpd}[$i] = $stuff[$i];}
#	if($on%100000==0){last;}
#$on++;
}



open(OUTPUT, ">", "MERGED_CPD_DB.txt")||die;
print OUTPUT "cpd\tsrc\tform\tmass\tchar\ttcdb\tname\tkegg\tchebi\thmdb\tpubch\tinchi\tbioc\n";
foreach $cpd (sort(keys %HoA)){
	print OUTPUT "$cpd\t"; 
	foreach $i (1 .. 12){
		print OUTPUT "$HoA{$cpd}[$i]\t";
	}
	print OUTPUT "\n";
}
	
