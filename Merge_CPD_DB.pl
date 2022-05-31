#use warnings;
use Sort::Naturally 'nsort';



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
qx{wget -N https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py};
open(INTRCH,    "getSubstrates.py")             ||die "unable to open getSubstrates.py: $!\n";


#GET TCDB - CHEBI DATA
$time=localtime;
print "INPUT TRANSPORTER SUBSTRATES $time\n";
while(<INTRCH>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\n\r]+//;
        (my $tcdb, my $chebs)=split("\t", $_);
        if($tcdb!~/TCDB/){$tcdb="TCDB:".$tcdb;}
        @CHEBS=();
        @CHEBS = ( $chebs =~ /(CHEBI.\d+)/g );
        foreach my $chebi (@CHEBS){ if($chebi=~/\d/ && $tcdb=~/\d/){$CHEB_TCDB{$chebi}{$tcdb}=1;}}
}



print "INPUT PATHBANK\n";
$on=0;
while(<INPB>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}
$time=localtime;
print "INPUT PUBCHEM on $on time $time\n";
while(<INPC>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}
$time=localtime;
print "INPUT HMDB on $on time $time\n";
while(<INHM>){	
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}
$time=localtime;
print "INPUT KEGG on $on time $time\n";
while(<INKG>){	
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}

$time=localtime;
print "INPUT BIOCYC on $on time $time\n";
while(<INBC>){	
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}
$time=localtime;
print "INPUT CHEBI on $on time $time\n";
while(<INCH>){	
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\n\r]+//;
	@stuff=split("\t", $_);
	$cpd=$stuff[0];		$CPDS{$cpd}=1;
	$src=$stuff[4];		$CPDSRC{$cpd}=$src;
	$form=$stuff[1]; 	$CPDFORM{$cpd}=$form;
	$mass=$stuff[2];	$CPDMASS{$cpd}=$mass;
	$char=$stuff[3];	$CPDCHAR{$cpd}=$char;
	if($on%100000==0){$time=localtime; print "on $on time $time cpd $cpd src $src form $form mass $mass char $char\n";}$on++;
	$tcdbs=$stuff[5];	@TCDBS=split(";",$tcdbs); foreach my $id (@TCDBS){ 	if($id=~/TCDB\:[\.\w]+/){	$CPD_TCDBS{$cpd}{$id}=1;}}
	$names=$stuff[6];	@NAMES=split(";",$names); foreach my $id (@NAMES){ 	if($id=~/\w{3}/){		$CPD_NAMES{$cpd}{$id}=1;}
										   		else{			print "cpd $cpd bad name $id\n";}}
	$keggs=$stuff[7];	@KEGGS=split(";",$keggs); foreach my $id (@KEGGS){ 	if($id=~/C\d+/){	 	$CPD_ALTS{$cpd}{$id}="K";}
												else{			print "cpd $cpd bad kegg $id\n";}}
	$chebs=$stuff[8];	@CHEBS=split(";",$chebs); foreach my $id (@CHEBS){ 	if($id=~/CHEBI\:\d+/){		$CPD_ALTS{$cpd}{$id}="C";}
												else{                   print "cpd $cpd bad cheb $id\n";}}
	$hmdbs=$stuff[9];	@HMDBS=split(";",$hmdbs); foreach my $id (@HMDBS){ 	if($id=~/HMDB\d+/){	 	$CPD_ALTS{$cpd}{$id}="H";}
												else{                   print "cpd $cpd bad hmdb $id\n";}}
	$pubcs=$stuff[10];	@PUBCS=split(";",$pubcs); foreach my $id (@PUBCS){ 	if($id=~/CID\:\d+/){	 	$CPD_ALTS{$cpd}{$id}="P";}
												else{                   print "cpd $cpd bad pubc $id\n";}}
	$inchs=$stuff[11];	@INCHS=split(";",$inchs); foreach my $id (@INCHS){ 	if($id=~/INCHI\:[A-Z\-]+/){	$CPD_ALTS{$cpd}{$id}="I";}
												else{                   print "cpd $cpd bad inch $id\n";}}
	$biocs=$stuff[12];	@BIOCS=split(";",$biocs); foreach my $id (@BIOCS){ 	if($id=~/\w/){			$CPD_ALTS{$cpd}{$id}="B";}
										   		else{			print "cpd $cpd bad bioc $id\n";}}
}

open(OUTPUT, ">", "MERGED_CPD_DB.txt")||die;
print OUTPUT "cpd\tsrc\tform\tmass\tchar\ttcdb\tname\tkegg\tchebi\thmdb\tpubch\tinchi\tbioc\n";
foreach $cpd (sort(keys %CPDS)){ 
	$src =$CPDSRC{$cpd};
	$form=$CPDFORM{$cpd};
	$mass=$CPDMASS{$cpd};
	$char=$CPDCHAR{$cpd};
	$kegg=''; $inchi=''; $pubch=''; $hmdb=''; $bioc=''; $chebi=''; $name='';
	%FORMS=(); %MASS=(); %CHARS=(); %NAMES=(); @NAMES=();
	#REMOVE BAD MATCH ALTS--FORM MASS CHAR NOT AGREE WITH MAIN CPD
	#ELSE IF NO FORM?MASS?CHAR FILL IN IF SAME VAL FOR MULTIPLE ALTS
        foreach my $alt (sort(keys %{$CPD_ALTS{$cpd}})){

		#formula
		if($form=~/\w/){ 
			if($CPDFORM{$alt} =~ /\w/ && $CPDFORM{$alt} ne $form){ delete($CPD_ALTS{$cpd}{$alt}); delete($CPD_ALTS{$alt}{$cpd}); next;}
		}
		else{ 	if($CPDFORM{$alt} =~ /\w/){ $FORMS{$CPDFORM{$alt}}++; }}
		#charge
		if($char=~/\w/){
		        if($CPDCHAR{$alt} =~ /\w/ && $CPDCHAR{$alt} ne $char){ delete($CPD_ALTS{$cpd}{$alt}); delete($CPD_ALTS{$alt}{$cpd}); next;}
		}
		else{   if($CPDCHAR{$alt} =~ /\w/){ $CHARS{$CPDCHAR{$alt}}++; }}
		#mass
		if($mass=~/\w/){
			#testing if masses within 0.1% mass - some DBs show diffent (to computer) but actually same mass
			if($CPDMASS{$alt} =~ /\w/){ #? should just use +/- 1, since 1 atomic mass equals proton?? 
				if($CPDMASS{$alt} > $mass+$mass*0.001 || $CPDMASS{$alt} < $mass-$mass*0.001){
					delete($CPD_ALTS{$cpd}{$alt}); delete($CPD_ALTS{$alt}{$cpd}); next;
		}	}	}
		else{   if($CPDMASS{$alt} =~ /\w/){ $MASS{$CPDMASS{$alt}}++; }}
		#names
		foreach my $name (keys %{$CPD_NAMES{$alt}}){ $NAMES{$name}++; }
	}


	if(keys %FORMS == 1){ foreach my $altf (keys %FORMS){ if($FORMS{$altf}>1 && $form!~/\w/){ $form=$altf; }}}
        if(keys %CHARS == 1){ foreach my $altf (keys %CHARS){ if($CHARS{$altf}>1 && $char!~/\w/){ $char=$altf; }}}
        if(keys %MASS == 1){  foreach my $altf (keys %MASS){  if($MASS{$altf}>1  && $mass!~/\w/){ $mass=$altf; }}}

	foreach my $name (keys %{$CPD_NAMES{$cpd}}){ $NAMES{$name}++; }
	foreach my $name (keys %NAMES){ push(@NAMES,$name); } 
	@NAMES=BestName(@NAMES);
	$name=join(";",@NAMES);

	#sum alts
	@KEGG=(); @PUBCM=(); @BIOCYC=(); @INCHI=(); @HMDB=(); @CHEBI=(); @TCDB=();
	foreach my $alt (sort(keys %{$CPD_ALTS{$cpd}})){
		if($CPD_ALTS{$cpd}{$alt} eq "C"){ 
			push(@CHEBI,$alt);
			foreach my $tcdb (sort(keys %{$CHEB_TCDB{$alt}})){push(@TCDB, $tcdb);}
		}
		if($CPD_ALTS{$cpd}{$alt} eq "H"){ push(@HMDB,$alt); }
		if($CPD_ALTS{$cpd}{$alt} eq "I"){ push(@INCHI,$alt); }
		if($CPD_ALTS{$cpd}{$alt} eq "B"){ push(@BIOCYC,$alt); }
		if($CPD_ALTS{$cpd}{$alt} eq "K"){ push(@KEGG,$alt); }
		if($CPD_ALTS{$cpd}{$alt} eq "P"){ push(@PUBCM,$alt); }
	}
	$kegg	=join(";",@KEGG);
	$inchi	=join(";",@INCHI);
	$pubch	=join(";",@PUBCM);
	$hmdb	=join(";",@HMDB);
	$bioc	=join(";",@BIOCYC);
	$chebi	=join(";",@CHEBI);
	$tcdb	=join(";",@TCDB);

	print OUTPUT "$cpd\t$src\t$form\t$mass\t$char\t$tcdb\t$name\t$kegg\t$chebi\t$hmdb\t$pubch\t$inchi\t$bioc\n";
}


#SUBROUTINES

sub CleanName{
        $nameX = $_[0];
        #remove junk punctuation/standardize
        $sta=0; $end=1;
        while($end ne $sta){
                $sta=$nameX;
                #swap greek symbols for text
                for my $g (0..$#GREEKL){        #fix pathbank and other greek symbols
                        if($nameX =~/($GREEKS[$g])/){
                                $nameX =~ s/$GREEKS[$g]/$GREEKL[$g]/g;
                }       }
                $nameX =~ s/\%2B(\d*)/$1\+/g;   #fix html +/- code (HMDB db)
                $nameX =~ s/\%2D(\d*)/$1\-/g;   #fix html +/- code (HMDB db)

                $nameX =~ s/\s+/_/g;
                $nameX =~ s/[^\w\-\+]+/_/g;
                $nameX =~ s/\_\+|\+\_/\+/g;
                $nameX =~ s/\_\-|\-\_/\-/g;
                $nameX =~ s/\-+/\-/g;
                $nameX =~ s/\++/\+/g;
                $nameX =~ s/\++\-+|\-+\++/\+/g;
                $nameX =~ s/\_+/\_/g;
                $nameX =~ s/(^[\_\W]+|[\_\W]+$)//g;

                #clear out junk descriptors
                $nameX =~ s/^(LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)//g;
                $nameX =~ s/^(UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)//g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)$//g;
                $nameX =~ s/[\b\_](UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)$//g;

                $end=$nameX;
        }
        return($nameX);
}

#REDUCE NUMBER OF ALT NAMES, PICK BEST
sub BestName{
        @NM = @_;
        #$kc=@NM;

        #GET LONGEST LETTER STRETCH, GET NON_WORD COUNTS
        %ODD=(); %LEN=();
        foreach my $name (@NM){
                $name=CleanName($name);
                $alt=$name;
                $alt=~s/[\W\_]+/\_/g;
                #get rid of runs of numbers
                while($alt=~/[\_\-]+(\d+[\-\_]\d+|\d{1,3}[A-Z]{1,2}\d{0,3}|\d{0,3}[A-Z]{1,2}\d{1,3})[\_\-]+/){
                        $alt=~s/[\_\-]+(\d+[\-\_]\d+|\d{1,3}[A-Z]{1,2}\d{0,3}|\d{0,3}[A-Z]{1,2}\d{1,3})[\_\-]+/\_/; }
                        #remove middle junk of numbers or numbers and letters
                while($alt=~/[\_\-]+(\d+|[A-Z]|\d{0,3}[A-Z]{1,2}\d{1,3}|\d{1,3}[A-Z]{1,2}\d{0,3})$/){
                        $alt=~s/[\_\-]+(\d+|[A-Z]|\d{1,3}[A-Z]{1,2}\d{0,3}|\d{0,3}[A-Z]{1,2}\d{1,3})$//; }
                        #remove junk after name ex: CARDIOLIPINS_20_4_20_1_18_2_18_2
                while($alt=~/^(\d+|[A-Z]|\d{1,3}[A-Z]{1,2}\d{0,3}|\d{0,3}[A-Z]{1,2}\d{1,3})[\_\-]+/){
                        $alt=~s/^(\d+|[A-Z]|\d{1,3}[A-Z]{1,2}\d{0,3}|\d{0,3}[A-Z]{1,2}\d{1,3})[\_\-]+//; }
                        #remove junk before name ex: 1_1_2_DI_5Z_8Z_11Z_14Z_EICOSATETRAENOYL_SN_GLYCERO_3_PHOSPHO_3_1_9Z_11Z_OCTADECADIENOYL_2_9Z_HEXADECENOYL_S>
                if($alt!~/[A-Z]{7,}/){$alt=$name; $alt=~s/[\W\_]+/\_/g;} #name too short after remove junk
                @PARTS = split("\_", $alt);
                $mp=0;
                foreach $p (@PARTS){
                        $p=~/([A-Z]+)/;
                        if(length($1)>$mp){ $mp=length($1); }
                }
                #$skip=0; foreach my $an (keys %LEN){if($an =~ /$alt/ && $an ne $alt){ $skip=1; last;}}
                #if($skip==1){next;}
                $LEN{$alt}=$mp; #whats longest letter stretch
                $mo=@PARTS;
                $ODD{$alt}=$mo; #how many name parts
        }

        #CHECK FOR NAMES >= 7 && <=30 CHARACTER STRINGS
        #there are some named with their protein/peptide == long
        #and there are some that are total junk: TG_A-17_0_I-18_0_I-15_0_RAC
        %GOOD=(); %SHORT=(); %LONG=();
        foreach my $name (keys %LEN){
                if($name !~/[A-Z]/){next;}
                   if($LEN{$name}>=7 && length($name) < 50){$GOOD{$name} =1; }
                elsif($LEN{$name}< 7 && length($name) < 50){$SHORT{$name}=1; }
                elsif(length($name) >= 50){                 $LONG{$name} =1; }
                else{ print "badname $name\n"; }
        }


        #COMBINE NESTED NAMES
        @NAMES=();
        foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %GOOD){
                $do=1;
                for my $i (0..$#NAMES){ #loop through added names, remove nested
                        if($NAMES[$i]=~/$name/){$do=0; last;}                   #the name is contained in another name, skip/redundant
                        if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
                }
                if($do==1){push(@NAMES,$name);}                                 #no name containment
                if($NAMES[9]=~/\w/){last;}                                      #limit to 10 good names
        }
        foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %LONG){
                $do=1;
                for my $i (0..$#NAMES){                                         #excluding nested names
                        if($NAMES[$i]=~/$name/){$do=0; last;}                   #the name is contained in another name, skip/redundant
                        if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
                }
                if($do==1){push(@NAMES,$name);}                                 #no name containment 
                if($NAMES[10]=~/\w/){last;}                                     #if 10 "good" names use only 1 long
        }
       foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %SHORT){
                $do=1;
                for my $i (0..$#NAMES){                                         #excluding nested names
                        if($NAMES[$i]=~/$name/){$do=0; last;}                   #the name is contained in another name, skip/redundant
                        if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
                }
                if($do==1){push(@NAMES,$name);}                                 #no name containment 
                if($NAMES[11]=~/\w/){last;}                                      #if 10 "good" and 1 long names use only 1 short
        }
        @NAMES=nsort(@NAMES);
        return(@NAMES);
}
