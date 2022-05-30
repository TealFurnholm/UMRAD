#use warnings;
use Sort::Naturally 'nsort';




#~~~~~~~~~~~~~~ PROCESS ~~~~~~~~~~~~~~~#
#THIS SCRIPT STARTS BY DOWNLOADING THE LATEST BIOLOGICAL COMPOUNDS AND REACTIONS FROM A VARIETY OF PUBLIC DATABASES
#SIX COMPOUND [CPD] DATABASES (INCHI, HMDB, KEGG, CHEBI, PUBCHEM, BIOCYC) CAN BE LINKED TO METABOLOMICS DATA
#THREE REACTION [RXN] DATABASES (BIOCYC, KEGG, RHEA) CAN BE LINKED TO THE REACTANT/PRODUCT COMPOUNDS AND THE TRANSPORT/ENZYME PROTEINS
#FOR EACH CPD AND RXN DB THE OUTPUTS ARE:
#THE CPD MASS, CHARGE, FORMULA, SOURCE (REPOSITORY), TRANSPORTER (if known), and CROSS-LINKED DB IDS ARE OUTPUT AS COMPOUND INFORMATION
#THE RXN DIRECTION, SOURCE (REPOSITORY), 3 CROSS-LINKED RXN DB IDS, EACH REACTANT AND PRODUCT COMPOUNDS FOR THE 3 RXN DBS, THEIR LOCATIONS (IN/OUTSIDE OF CELL), AND DIRECTION OF CPD TRANSPORT.
#~~~~~
#PREVIOUS VERSION I TRIED TO USE THE VARIOUS DATABASES TOGETHER, TO FILL IN MISSING INFORMATION
#KEEPING THEM SEPARATE NOW BECAUSE TOO MANY ERRORS IN PUBLIC REPOSITORIES THAT AMPLIFY UPON COMBINATION
#YOU'LL SEE A LOT OF PARSING AND "BLANK" CHECKING CODE
#HOPEFULLY THE OUTPUTS FROM EACH REPOSITORY WILL BE RELIABLY COMBINABLE
#THE CPD AND RXN DBs WILL BE USED IN THE CREATION OF THE PROTEIN ALIGNMENT DATABASE FOR UMRAD: https://github.com/TealFurnholm/UMRAD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#




################################################
########    SET INPUT/OUTPUT FILES    ##########
################################################

#LIVE DOWNLOADS: -- leaving here so Anders can split how he needs in his snakemake while I can put this into GitHub for non-UMich users
$time=localtime;
print "DOWNLOAD INPUTS time $time\n";
qx{wget -N https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py};
$ec_up = `wget -q -O - https://ftp.expasy.org/databases/enzyme/enzyme.dat`;



#OPEN/CHECK INPUTS:
$time=localtime;
print "OPEN / CHECK INPUTS time $time\n";
open(INTRCH,	"getSubstrates.py")		||die "unable to open getSubstrates.py: $!\n";


#1. SET BIOCYC DIRECTORY
	$bioc_dir='/geomicro/data2/kiledal/UMRAD/Universal_Biological_Compounds_Database/BIOCYC_NF/*';
	#####$mioc_dir='/geomicro/data2/tealfurn/URDB/UNIPROT/FUNCTIONS/BIOCYC/*'; #!!!!##

#2. get metacyc - need to change script to load it first and not change if "exists",   
	# so all other records will be beholden to curated Metacyc data 
	# there were problems with LTR/RTL annotations in individual oranism files probably others
	# caused left/right cpds end up on both sides of rxn

	$mdir = $bioc_dir;
	#$mdir = $mioc_dir; #!!!!
	$mdir =~ s/\*$/meta/; 	#getting metacyc files first, eventually tweak so uses these first, 
				#not replaced by other less well annotated oranism files
	$cdatf=$mdir."/compounds.dat";
	$rdatf=$mdir."/reactions.dat";
	$mxmlf=$mdir."/metabolic-reactions.xml";
	$pseqf=$mdir."/protein-seq-ids.dat";

	#reaction-links files have rxn -> EC#: maybe later check how much improved annotataion if use this file too
	#RXN0-2522	EC-2.7.1.195
	#2.3.1.126-RXN	EC-2.3.1.126
	#PYRAZIN-RXN	EC-3.5.1.4

	if(-f $cdatf && -s $cdatf && -f $rdatf && -s $rdatf && -f $mxmlf && -s $mxmlf){ 
		#print "mdir $mdir\ncdatf $cdatf \nrdatf $rdatf \nmxmlf $mxmlf \nzdatf $zdatf \npdatf $pdatf\n";
		push(@CPD_DAT,$cdatf);
		push(@RXN_DAT,$rdatf);
		push(@METRXN_XML,$mxmlf);
		if(-s $pseqf){push(@PRT_SEQ,$pseqf);}
		else{	$pseqf=$mdir."/uniprot-seq-ids.dat";
			push(@PRT_SEQ,$pseqf);}
	}

#3. next get the rest of the organism directories
$x=''; 
$x = qx{ls -d $bioc_dir}; 
@BDIRS =split("\n", $x);
foreach my $bdir (@BDIRS){
	if($bdir =~ /\/meta$|\/meta\/$/i){ next;} #skip meta, already have added it
	$cdatf=$bdir."/compounds.dat";
	$rdatf=$bdir."/reactions.dat";
	$mxmlf=$bdir."/metabolic-reactions.xml";
	$pseqf=$bdir."/protein-seq-ids.dat";
	if(-f $cdatf && -s $cdatf && -f $rdatf && -s $rdatf && -f $mxmlf && -s $mxmlf){ 
		push(@CPD_DAT,$cdatf);
		push(@RXN_DAT,$rdatf);
		push(@METRXN_XML,$mxmlf);
		if(-f $pseqf && -s $pseqf){push(@PRT_SEQ,$pseqf);}
	}
	#$kc=@CPD_DAT;
	#if($kc>=10){last;} #!!!!
}
if($CPD_DAT[0]	 !~/\w/){ print "missing BIOCYC compounds.dat files @CPD_DAT.\n"; 		die; }
if($RXN_DAT[0]	 !~/\w/){ print "missing BIOCYC reactions.dat files @RXN_DAT.\n"; 		die; }
if($METRXN_XML[0]!~/\w/){ print "missing BIOCYC metabolic-reactions.xml files @METRXN_XML.\n";	die; }
$cd=@CPD_DAT;
$rd=@RXN_DAT;
$mx=@METRXN_XML;
$ps=@PRT_SEQ;
print "cd $cd rd $rd mx $mx ps $ps\n";

#OPEN/CHECK OUTPUTS:
open(OUTBCC, ">", "BIOCYC_CPD_DB.txt")||die;
#open(OUTBCC, ">", "BIOCYC_CPD_DB_test.txt")||die;
open(OUTBCR, ">", "BIOCYC_RXN_DB.txt")||die;
#open(OUTBCR, ">", "BIOCYC_RXN_DB_test.txt")||die; #!!!!
################################################
################################################








##################################################################################
######################   	LOAD UNIVERSAL INFO	 #########################
##################################################################################

#elements sorted long to short to remove in correct order
#@ELEMENTS = ("AC","AG","AL","AM","AR","AS","AT","AU","TM","BA","BE","BH","BI","BK","BR","TS","CA","CD","CE","CF","CL","CM","CN","CO","CR","CS","CU","DB","DS","DY","ER","ES","EU","XE","FE","FL","FM","FR","GA","GD","GE","YB","HE","HF","HG","HO","HS","ZN","IN","IR","ZR","KR","LA","LI","LR","LU","LV","MC","MD","MG","MN","MO","MT","NA","NB","ND","NE","NH","NI","NO","NP","OG","OS","PA","PB","PD","PM","PO","PR","PT","PU","RA","RB","RE","RF","RG","RH","RN","RU","SB","SC","SE","SG","SI","SM","SN","SR","TA","TB","TC","TE","TH","TI","TL","B","C","F","H","I","K","U","V","W","P","Y","N","S","O");
@GREEKS = ("α","β","γ","δ","ε","ζ","η","θ","ι","κ","λ","μ","ν","ξ","ο","π","ρ","ς","σ","τ","υ","φ","χ","ψ","ω");
@GREEKL = ("ALPHA","BETA","GAMMA","DELTA","EPSILON","ZETA","ETA","THETA","IOTA","KAPPA","LAMBDA","MU","NU","XI","OMICRON","PI","RHO","SIGMA","SIGMA","TAU","UPSILON","PHI","CHI","PSI","OMEGA");

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
	$TCDB{$tcdb}=1;
}
$tckc=keys %TCDB;
$tcdbkc=keys %CHEB_TCDB;
print "there are $tcdbkc chebis and tcdb $tckc\n";
undef(%TCDB);
#  ! consider aligning tcdb-chebi uniprot prots to uniprot db before now,
#  get the annotations to put here not just chebi but other cpd IDs...


# GET EC -> UNIPROT
$time=localtime;
print "INPUT ECs TO UPIDs $time\n";
@EC_UP = split("\n", $ec_up);
$inrec=0;
%UPID=();
foreach my $x (@EC_UP){
	@stuff=split("\t",$x);
	$x=uc($x);
	if($x=~/^ID\s+(\d+\.\d+\.\d+\.N*\d+)/){ $ec=$1; $inrec=1; $alt='';}
	#get single obsoletes - too frequently used to ignore
	#ignore obsoletes split into multiple ECs
	if($inrec==1 && $x=~/TRANSFERRED.ENTRY\:\s+([\w\.]+)\.$/){ $alt=$1; $EC2EC{$ec}=$alt; $EC2EC{$alt}=$ec;}
	if($inrec==1 && $x=~/^\s*DR\s/){
	        @fnd= ($x =~ /\s+(\w+)\,/g);
	        foreach my $upid (@fnd){
	                if($upid=~/\w/ && $ec=~/\w/){
	                        $UPID{$upid}=1;
	                        $EC2UPID{$ec}{$upid}=1;
	}       }       }
	if($x=~/^\s*\/\/\s*$/){ $inrec=0; }
}
#distribute upids to obsoletes
foreach my $ec (keys %EC2UPID){
	$alt = $EC2EC{$ec};
	if($alt !~/\w/){next;}
	foreach my $upid (keys %{$EC2UPID{$ec}}){ $EC2UPID{$alt}{$upid}=1;}
}
undef(%EC2EC);
##################################################################################
######################      DONE UNIVERSAL INFO          #########################
##################################################################################








#~~~~~~~~~~~~~~~~~~~~~~~~ LOAD REPOSITORY DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~#






##################################################################################
####################    BIOCYC COMPOUNDS AND REACTIONS      ######################
##################################################################################

#transport and rxn left/right in/out
#you were adding INSIDE/OUTSIDE to specific cpd metaids in the metabolic-reactions.xml file, both from cpd and rxn
#you can use cpd metaid or "compartment" store in $CMPD_LOC{$metacpd}=$loc;
#then use rxn product/reactant cpd metaids to get more inside/outside store in $CMPD_LOC{$metacpd}=$loc;
#then use the cpd loc and metacpd to store in/out in left/right cpds $LCPD{$CPD_METAID{$left}}=$CMPD_LOC{$left};
#then use %RCPD/LCPD to populate rxn left/right rxn cpds and in/out $LRXN_CPD{$rxn}{$cpd}=$LCPD{$cpd};
#finally use reactions.dat to get more left/right-in/out info
#then for each reaction determine if out/in differs for same cpd on left and right -> a transporter and the direction


#######################################################	
#FIRST LOAD ALL BIOCYC CPDS
#######################################################	
$time=localtime;
print "INPUT BIOCYC compounds.dat time $time\n";
$on=1;
$rec=@CPD_DAT;
foreach my $file (@CPD_DAT){
#print "on $on of $rec file $file time $time\n"; #!!!!

	if($on%1000==0){$time=localtime; print "on $on of $rec file $file time $time\n"; } $on++;
	open(INBIOC,$file)||die "unable to open $file: $!\n";
	$inrxn=0;
	while(<INBIOC>){
	        $_ = uc($_);
	        $_=~s/[\r\n]+//;
	        if($_ =~ /UNIQUE\-ID\s+\-\s+(\S+)/){ 
			%ALTS=(); $cpd=$1; $ALTS{$cpd}=1; %NAMES=();  %CHFO=(); 
			$at=''; $nu=''; $form=''; $mass=''; $name=''; $inrxn=1; next;}
		if($inrxn==1){
	                if($_ =~ /LIGAND\-CPD\s+\"(C\d+)\"/){		$alt=$1; 		$ALTS{$alt}=1;	next;}
			if($_ =~ /INCHIKEY\=([\w\-]+)/){		$alt="INCHI:".$1;	$ALTS{$alt}=1;	next;}
			if($_ =~ /DBLINKS\s+\-\s+\(CHEBI\s+\"(\d+)/){	$alt="CHEBI:".$1;	$ALTS{$alt}=1;	next;}
			if($_ =~ /DBLINKS\s+\-\s+\(PUBCHEM\s+\"(\d+)/){	$alt="CID:".$1;		$ALTS{$alt}=1;	next;}
	  		if($_ =~ /DBLINKS\s+\-\s+\(HMDB\s+\"(HMDB\d+)/){$alt=$1;		 
				while(length($alt)<11){$alt=~s/HMDB/HMDB0/;}			$ALTS{$alt}=1;	next;}
			if($_ =~ /MOLECULAR.WEIGHT\D+([\d\.]+)/){				$mass=$1;	next;}
			if($_ =~ /CHEMICAL.FORMULA\W+([A-Z]{1,2})\s+(\d+)/i){ $at=$1; $nu=$2;	$CHFO{$at}=$nu;	next;}
			#get name
	                if($_ =~ /TYPES\s+\-\s+(.*)/ && $name eq ''){	$name=$1;		$NAMES{$name}=1; next;}
	                if($_ =~ /COMMON\-NAME\s+\-\s+(.*)/){		$name=$1;		$NAMES{$name}=1; next;}
			if($_ =~ /SYNONYMS\s+\-\s+(\S.*)/){		$name=$1;		$NAMES{$name}=1; next;}
	        }
	        if($_=~/^\/\/\s*$/){
			$inrxn=0;
			if($cpd !~/\w/){next;}
			if($mass=~/\d/){$CMPD_MASS{$cpd}=$mass;}
			foreach my $name (keys %NAMES){
				#remove html code from name
				if($name=~/\<.*\>/){ $name=~s/\<[^\>]*\>//g; }
				if($name=~/\&[^\;]*\;/){
					foreach my $gl (@GREEKL){$name =~ s/\&$gl\;/$gl/ig;}
					$name=~s/\&\w*DASH\;/\-/ig;
					$name=~s/\&[^\;]*\;/\_/ig;
				}
				$name=CleanName($name);
				if($name=~/\w/ && $name!~/^(COMPOUNDS|ALL-COAS|RINGS|PSEUDO-COMPOUNDS|IONS|BASES)$/){
					$CMPD_NAME{$cpd}{$name}="BC";
				}
			}
			foreach my $alt (keys %ALTS){   
				if($alt=~/\w/){$CMPD_ALTS{$cpd}{$alt}="BC";}}
			$form='';
			foreach my $atom (sort(keys %CHFO)){ 
				$form.=$atom; 
				if($CHFO{$atom}!=1){$form.=$CHFO{$atom};}
			}
			$CMPD_FORM{$cpd}=CleanName($form);
	}	}	
}
#######################################################	
#######################################################	



#######################################################
#LOAD BICYC UPIDS
#######################################################
$on=1;
$rec=@PRT_SEQ;
print "INPUT BIOCYC prt-seq\n";
foreach my $file (@PRT_SEQ){
	$/="(";
	if($on%1000==0){$time=localtime; print "on prtseq $on of $rec file $file time $time\n"; } $on++;
	@UPIDS=();
	open(INPSQ, $file)||die;
	while(<INPSQ>){
		if($_!~/UNIPROT\:\w+/){next;}
		$_=uc($_);
		$_=~s/[\r\n]+//;
		@UPIDS = ( $_ =~ /UNIPROT:([^\"]+)/g );
		@stuff=split('\s+',$_);
		$rxn=$stuff[0];
		if($stuff[1]=~/(\d+\.\d+\.\d+\.\d+)/){ $ec=$1; } 
		else{$ec='';}
		foreach my $upid (@UPIDS){ 	$RXN_UPID{$rxn}{$upid}="BC";
						$UPID_RXN{$upid}{$rxn}="BC";
			if($ec=~/\d/){		$EC2UPID{$ec}{$upid}=1;
						$RXN_EC{$rxn}{$ec}="BC";
			}
	}	}
	$/="\n"; #reset delimiter
}
$upkc=keys %UPID_RXN;
$rxkc=keys %RXN_UPID;
print "upkc $upkc rxkc $rxkc\n";
#######################################################	
#######################################################	







#######################################################	
#NOW LOAD BIOCYC CPD-RXN XML FILES
#######################################################	
$time=localtime;
print "INPUT BIOCYC metabolic-reactions.xml time $time\n";
$on=1;
$rec=@METRXN_XML;
foreach my $file (@METRXN_XML){
#	print "on $on of $rec file $file time $time\n"; #!!!!

	$on++;
	if($on%1000==0){$time=localtime; print "on $on of $rec file $file time $time\n"; }
        open(INXML,$file)||die "unable to open $file: $!\n";
	$incmpd=0;
	BIOCXML: while(<INXML>){
	        if($_ !~ /\w/){next;}
	        $_ = uc($_);
		$_ =~ s/^[\n\r\s]+|[\n\r\s]+$//g;
		if($_ =~ /^\</){ $line=$_; if( $line !~ /\>$/ ){next BIOCXML;} }
		while($line !~ /\>$/){ $line.=" ".$_; 		
			if( $line !~ /\>$/ ){next BIOCXML;}
			else{last;}
		}
		$line=~s/\s+//g;

	       	#COMPOUND INFO
	       	if($line=~/SPECIESMETAID\=\"([^\"]+)/i){
			$incmpd=1; $metacpd=$1; $cpd=''; %ALTS=(); $name=''; $form=''; $char=''; $ch='';
			$loc="INSIDE"; if($metacpd =~ /\_(OUT|E|P)$/i){$loc="OUTSIDE";}
		}
	       	if($incmpd==1){
			if($line=~/name\=\"([^\"]+)/i){				$name=CleanName($1);}
			if($line=~/compartment\=\"([^\"]+)/i){			$comp=$1; 
				if($comp =~ /^(CCO[\d\_]+OUT|E|P)$/i){		$loc="OUTSIDE";}}
	       	        if($line=~/CHEMICALFORMULA[^\"]+\"([^\"]+)/i){ 		$form=CleanName($1);}
	       	        if($line=~/fbc.charge\=\"([\-]*)(\d+)/i){ 		$ch=$2.$1;
				#FBC:CHARGE="0"
				if($ch !~ /\-/ && $ch !~ /^0$/){$ch.="+";} 	$char=$ch;}
	       	        if($line=~/identifiers.org.biocyc.[^\:]+\:([^\"]+)/i){	$cpd=$1; 		$ALTS{$1}=1;}
	       	        if($line=~/identifiers.org.kegg.compound.(C\d+)/i){ 				$ALTS{$1}=1;}	
	       	        if($line=~/identifiers.org.chebi.(CHEBI.\d+)/i){				$ALTS{$1}=1;}
			if($line=~/identifiers.org.inchikey.([A-Z\-]+)/i){    	$id="INCHI:".$1;	$ALTS{$id}=1;}
			if($line=~/identifiers.org.pubchem.compound.(\d+)/i){ 	$id="CID:".$1;		$ALTS{$id}=1;}
			if($line=~/identifiers.org.hmdb.(HMDB\d+)/i){ 		$id=$1; 
				while(length($id)<11){$id=~s/HMDB/HMDB0/;}				$ALTS{$id}=1;}
		}
	       	if($line=~/\<\/SPECIES\>/i){
	       	        $incmpd=0;
	       	        if($cpd !~/\w/){next;}
			$CPD_METAID{$metacpd}=$cpd;
			$CMPD_LOC{$metacpd}=$loc;
			$CMPD_NAME{$cpd}{$name}="BC";
			$CMPD_FORM{$cpd}=$form;
			$CMPD_CHAR{$cpd}=$char;
			$CMPD_ALTS{$cpd}{$cpd}="BC";
			#print "cpd $cpd loc $loc char $char form $form meta $metacpd name $name\n";
			foreach my $alt (keys %ALTS){if($alt =~/\w/){ $CMPD_ALTS{$cpd}{$alt}="BC";}}
		}
		#RESET LINE
		$line='';
	} 


	#NOW GET RXN INFO - MUST BE DONE IN THIS ORDER AND SEPARATELY FROM COMPOUNDS FROM THIS XML FILE
	if($on%1000==0){$time=localtime; print "on $on of $rec file $file time $time\n"; }
	open(INXML, $file)||die "unable to open $file: $!\n";
	$inrxn=0;
	BIOCRXML: while(<INXML>){
	        if($_ !~ /\w/){next;}
	        $_ = uc($_);
		$_ =~ s/^[\n\r\s]+|[\n\r\s]+$//g; #DELETE LEADING/TRAILING WHITESPACE
		if($_ =~ /^\</){ $line=$_; if( $line !~ /\>$/ ){next BIOCRXML;}	}
		while($line !~ /\>$/){ 	$line.=" ".$_; 		
			if( $line !~ /\>$/ ){next BIOCRXML;}
			else{last;}
		}
		$line=~s/\s+//g;


	
		#REACTION INFO
		if($line=~/reactionmetaid\=\"([^\"]+)/i){
			$inrxn=1; 
			$metarxn=$1; 
			$rxn=''; $ec='';
			$dir="LTR"; 
			%LCPD=(); %RCPD=(); %RXNIDS=();
		}
		if($line=~/\<LISTOFREACTANTS\>/i){	$reacts=1; $prods=0;}
		if($line=~/\<listOfProducts/i){ 	$reacts=0; $prods=1;}
		if($line=~/\<\/LISTOFREACTANTS\>/i){	$reacts=0;	    }
		if($line=~/\<\/listOfProducts/i){ 	$prods=0;	    }
		if($inrxn==1){
			if($line=~/identifiers.org.biocyc.[^\:]+\:([^\"]+)/i){		$rxn=$1;	$RXNIDS{$rxn}=1;}
			if($line=~/identifiers.org.rhea.(\d+)/i){			$rhea=$1; 	$RXNIDS{$rhea}=1;}
			if($line=~/identifiers.org.kegg.reaction.(R\d+)/i){		$kegg=$1; 	$RXNIDS{$kegg}=1;}
			if($line=~/REVERSIBLE\W+true/i){				$dir="BOTH";	}
			if($line=~/ers.org.ec.code.(\d+\.\d+\.\d+\.\d+[\.\d]*)/i){	$ec=$1;	  	}
			if($reacts==1 && $line=~/SPECIES\=\"([^\"]+)/i){		
				$left=$1;
				if($CMPD_LOC{$left}!~/\w/){
					if($left=~/\_(OUT|E|P)$/){$CMPD_LOC{$left}="OUTSIDE";} 
					else{$CMPD_LOC{$left}="INSIDE";}}	
				if($CPD_METAID{$left} =~/\w/){
					if($CMPD_LOC{$left}=~/\w/){$LCPD{$CPD_METAID{$left}}=$CMPD_LOC{$left};}
					else{$LCPD{$CPD_METAID{$left}}="INSIDE";}}
			}
			if($prods==1 && $line=~/SPECIES\=\"([^\"]+)/i){
		 		$right=$1;
				if($CMPD_LOC{$right}!~/\w/){
					if($right=~/\_(OUT|E|P)$/){$CMPD_LOC{$right}="OUTSIDE";} 
					else{$CMPD_LOC{$right}="INSIDE";}}	
				if($CPD_METAID{$right}=~/\w/){
					if($CMPD_LOC{$right}=~/\w/){$RCPD{$CPD_METAID{$right}}=$CMPD_LOC{$right};}
					else{$RCPD{$CPD_METAID{$right}}="INSIDE";}}
		}	}
		if($_=~/\<\/REACTION\>/i){
			$inrxn=0;
			if($rxn !~ /\w/ || $metarxn !~ /\w/){next;}
			if($ec=~/\d/){$RXN_EC{$rxn}{$ec}="BC";}
			$RXN_DIR{$rxn}=$dir;
			foreach my $cpd (keys %LCPD){ $LRXN_CPD{$rxn}{$cpd}=$LCPD{$cpd};	$LCPD_RXN{$cpd}{$rxn}=$LCPD{$cpd};}
			foreach my $cpd (keys %RCPD){ $RRXN_CPD{$rxn}{$cpd}=$RCPD{$cpd};	$RCPD_RXN{$cpd}{$rxn}=$RCPD{$cpd};}
			foreach my $alt (keys %RXNIDS){ if($alt=~/\w/){	$RXN_ALTS{$rxn}{$alt}="BC";}}			
		}
	
		#RESET LINE
		$line='';
	}
	undef(%CPD_METAID);
	undef(%CMPD_LOC);
}
#######################################################	
#######################################################	






#######################################################	
#NOW LOAD BIOCYC REACTION DATA FILES
#######################################################	
$time=localtime;
print "INPUT BIOCYC reactions.dat time $time\n";
$on=1;
$rec=@RXN_DAT;
foreach my $file (@RXN_DAT){
	#print "on rxdat $on of $rec file $file time $time\n"; #!!!!

	if($on%1000==0){$time=localtime; print "on $on of $rec file $file time $time\n"; } $on++;
	open(RXDAT, $file)||die "unable to open $file:$!\n";
	$inrxn=0;
	while(<RXDAT>){
		$_ = uc($_);
		$_=~s/[\r\n]+//;
		if($_ =~ /UNIQUE-ID\s+\-\s+(\S+)/){$rxn=$1; %RXNIDS=(); $RXNIDS{$rxn}=1; $inrxn=1; 
			$dir=''; $ec=''; $up=''; $rhea=''; $out=0; $in=0; $kegg=''; %LEFT=(); %RIGHT=(); %UP=(); next;}
		if($inrxn==1){
			#GET L/R THEN CHECK IF THERE IS AN INDICATION OF INSIDE OR OUTSIDE OF CELL
			#L/R ARE ACTUAL CPD IDs not metaids
			if($_ =~ /LEFT\s\-\s(\S+)/){	$inleft=1;  $left=$1;	$LEFT{$left}="INSIDE"; 	next;}
			if($inleft==1){ 
				if($_ =~ /\^COMPARTMENT\s\-\s(CCO-MIDDLE|CCO-OUT)$/){$LEFT{$left}="OUTSIDE"; $out=1;}
				$inleft=0;  $left='';
			}
			if($_ =~ /RIGHT\s\-\s(\S+)/){	$inright=1; $right=$1;	$RIGHT{$right}="INSIDE"; next;}
			if($inright==1){
				if($_ =~ /\^COMPARTMENT\s\-\s(CCO-MIDDLE|CCO-OUT)$/){$RIGHT{$right}="OUTSIDE"; $out=1;}
				$inright=0; $right='';
			}
	
			#GET REMAINING INFO
			if($_ =~ /DBLINKS.*RHEA\s+\"(\d+)\"/){          $rhea=$1; $RXNIDS{$rhea}=1; 	next;}
			if($_ =~ /LIGAND\-RXN\s+\"(R\d+)\"/){           $kegg=$1; $RXNIDS{$kegg}=1; 	next;}
			if($_ =~ /EC\-NUMBER\s\-\sEC\-(\d+\.\d+\.\d+\.\d+[\.\d]*)\s*$/){$ec=$1;	 	next;}
			if($_ =~ /UNIPROT\s+\"([^\"]+)\"/){		$up=$1;	  $UP{$up}=1;		next;}
			if($_ =~ /REACTION.DIRECTION.*(LEFT.*RIGHT)/){	$dir="LTR";  			next;}
			if($_ =~ /REACTION.DIRECTION.*(RIGHT.*LEFT)/){	$dir="RTL";  			next;}
			if($_ =~ /REACTION.DIRECTION.*(REVERSIBLE)/){	$dir="BOTH";  			next;}
	
			#CHECK IF ANY INSIDE COMPOUNDS
			if($_=~ /RXN.LOCATIONS\s\-\s(\S+)/){
				$loc=$1;
				   if($loc =~ /^(CCO-PERIPLASM|CCO-PERI-BAC-CCO-PERI-BAC|CCO-PERI-BAC|CCO-OUT-CCO-EXTRACELLULAR|CCO-OUT)$/){$out=1;}
				elsif($loc =~ /^(CCO-EXTRACELLULAR-CCO-UNKNOWN-SPACE|CCO-EXTRACELLULAR|CCI-PERI-BAC-GN)$/){$out=1;}
				else{$in=1;} 
				next;
			}
		}
		if($_=~/^\/\/\s*$/){
			$inrxn=0;
			if($rxn!~/\w/){next;}
			foreach my $alt (keys %RXNIDS){ if($alt=~/\w/){ $RXN_ALTS{$rxn}{$alt}="BC";}}
	
			#FIRST FLIP RTL RXN CPDS SO ALL ARE LTR
			if($dir eq "RTL"){
				%RSWITCH=(); %LSWITCH=();
				foreach my $left (keys %LEFT){ 	$LSWITCH{$left}=$LEFT{$left};}
				foreach my $right (keys %RIGHT){$RSWITCH{$right}=$RIGHT{$right};}
				undef(%RIGHT); undef(%LEFT);
				foreach my $xcpd (keys %LSWITCH){$RIGHT{$xcpd}=$LSWITCH{$xcpd};}
				foreach my $xcpd (keys %RSWITCH){$LEFT{$xcpd} =$RSWITCH{$xcpd};}
				$dir="LTR";
			}
	
			#NOW CHECK L/R IF EXISTS IN L/R_CMPD ALREADY, MORE RELIABLE
			#OR IF RXN HAS ONLY ONE RXN-LOCATION
			foreach my $cpd (keys %LEFT){
				if($LRXN_CPD{$rxn}{$cpd}=~/INSIDE|OUTSIDE/){ $LEFT{$cpd}=$LRXN_CPD{$rxn}{$cpd}; } #USE xml data, more specific
				elsif($in==1 && $out==0){$LEFT{$cpd}="INSIDE";} #all cpds are in
				elsif($in==0 && $out==1){$LEFT{$cpd}="OUTSIDE";} #all cpds are out
				else{}
				$LRXN_CPD{$rxn}{$cpd}=$LEFT{$cpd};
				$LCPD_RXN{$cpd}{$rxn}=$LEFT{$cpd};
			}
			foreach my $cpd (keys %RIGHT){
				if($RRXN_CPD{$rxn}{$cpd}=~/INSIDE|OUTSIDE/){ $RIGHT{$cpd}=$RRXN_CPD{$rxn}{$cpd}; } #USE xml data, more specific
				elsif($in==1 && $out==0){$RIGHT{$cpd}="INSIDE";} #all cpds are in
				elsif($in==0 && $out==1){$RIGHT{$cpd}="OUTSIDE";} #all cpds are out
				else{}
				$RRXN_CPD{$rxn}{$cpd}=$RIGHT{$cpd};
				$RCPD_RXN{$cpd}{$rxn}=$RIGHT{$cpd};
			}
	
			#OTHER INFO
			if($dir =~ /(LTR|BOTH)/){$RXN_DIR{$rxn}=$dir;}
					    else{$RXN_DIR{$rxn}="LTR";}
			if($ec=~/\d/){ $RXN_EC{$rxn}{$ec}="BC"; }
			foreach my $upid (keys %UP){
				if($upid=~/\w/){$UPID_RXN{$upid}{$rxn}="BC"; 
						$RXN_UPID{$rxn}{$upid}="BC";
			}	}
	
			#check for transport
			foreach my $cpd (keys %{$LRXN_CPD{$rxn}}){ #check whether inside->outside or viceversa
				#cpd must be on both sides out/in, so only need to check one rxn side
				if(exists( $RRXN_CPD{$rxn}{$cpd}) && $RRXN_CPD{$rxn}{$cpd} ne $LRXN_CPD{$rxn}{$cpd}){
					if($RRXN_CPD{$rxn}{$cpd} eq "INSIDE"  && $LRXN_CPD{$rxn}{$cpd} eq "OUTSIDE"){
					    	$TRANS_RXN{$rxn}{$cpd}="IMPORT"; 
					    	$TRANS_CPD{$cpd}{$rxn}="IMPORT";}
					if($RRXN_CPD{$rxn}{$cpd} eq "OUTSIDE" && $LRXN_CPD{$rxn}{$cpd} eq "INSIDE"){ 
					    	$TRANS_RXN{$rxn}{$cpd}="EXPORT"; 
					    	$TRANS_CPD{$cpd}{$rxn}="EXPORT";}
					if($dir eq "BOTH"){
						$TRANS_RXN{$rxn}{$cpd}="BIPORT"; 
						$TRANS_CPD{$cpd}{$rxn}="BIPORT";}
}	}	}	}	}
#######################################################	
#######################################################	




#######################################################	
#OUTPUT BIOCYC CPDs
#######################################################	
$time=localtime;
print "OUTPUT BioCYC CPDS time $time\n";
print OUTBCC "cpd\tformula\tmass\tcharge\tdb_src\ttcdbs\tnames\tkeggcpd\tschebcpd\thmdbcpd\tpubccpd\tinchcpd\tbioccpd\n";
foreach my $cpd (sort(keys %CMPD_ALTS)){
	$form=''; $char=''; $mass=''; $name='';
	$keggcpd='';   $chebcpd='';   $hmdbcpd='';
	$pubccpd='';   $inchcpd='';   $bioccpd='';
	$form=$CMPD_FORM{$cpd};
	$char=$CMPD_CHAR{$cpd};
	$mass=$CMPD_MASS{$cpd};
	$sorc=$CMPD_ALTS{$cpd}{$cpd};

	#FIX NAMES
	@NAMES=();
	foreach my $id (keys %{$CMPD_NAME{$cpd}}){push(@NAMES, $id);}
	@NAMES=BestName(@NAMES);
	$name=join(";",@NAMES);

	#sort identifiers by database CHECK EVERY ALT CPD FORMULA MATCHES
	foreach my $alt (sort(keys %{$CMPD_ALTS{$cpd}})){
		$alf=$CMPD_FORM{$alt};
		if($alf ne $form && $alf=~/\w/){ #two formulas dont match, remove alt
			   if($alt=~/^C\d+$/){		$keggcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			elsif($alt=~/^CHEBI.\d+$/){	$chebcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			elsif($alt=~/^HMDB\d+$/){	$hmdbcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			elsif($alt=~/^CID.\d+$/){	$pubccpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			elsif($alt=~/^INCHI.\S+$/){	$inchcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			else{				$bioccpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
			delete($CMPD_ALTS{$cpd}{$alt}); 
			next;
		}
		   if($alt=~/^C\d+$/){		if($keggcpd!~/$alt/){$keggcpd.=$alt.";"; }}
		elsif($alt=~/^CHEBI.\d+$/){	if($chebcpd!~/$alt/){$chebcpd.=$alt.";"; }}
		elsif($alt=~/^HMDB\d+$/){	if($hmdbcpd!~/$alt/){$hmdbcpd.=$alt.";"; }}
		elsif($alt=~/^CID.\d+$/){	if($pubccpd!~/$alt/){$pubccpd.=$alt.";"; }}
		elsif($alt=~/^INCHI.\S+$/){	if($inchcpd!~/$alt/){$inchcpd.=$alt.";"; }}
		else{				if($bioccpd!~/\Q$alt\E/){$bioccpd.=$alt.";"; }}
	}
	$keggcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $keggcpd); @CPDS=nsort(@CPDS); $keggcpd=join(";",@CPDS);
	$chebcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $chebcpd); @CPDS=nsort(@CPDS); $chebcpd=join(";",@CPDS);
	$cc=@CPDS; %TCDBS=(); @TCDBS=();  $tcdb=''; #GET CHEBI TCDBS
	if($cc<=3){foreach my $ch (@CPDS){ foreach my $tc (sort(keys %{$CHEB_TCDB{$ch}})){ $TCDBS{$tc}=1; }}
	if(keys %TCDBS > 0){ 
		foreach my $tc (sort(keys %TCDBS)){ push(@TCDBS,$tc); } @TCDBS=nsort(@TCDBS); $tcdb=join(";",@TCDBS);}}
	$hmdbcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $hmdbcpd); @CPDS=nsort(@CPDS); $hmdbcpd=join(";",@CPDS);
        $pubccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $pubccpd); @CPDS=nsort(@CPDS); $pubccpd=join(";",@CPDS);
	$inchcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $inchcpd); @CPDS=nsort(@CPDS); $inchcpd=join(";",@CPDS);
	$bioccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $bioccpd); @CPDS=nsort(@CPDS); $bioccpd=join(";",@CPDS);
	print OUTBCC "$cpd\t$form\t$mass\t$char\t$sorc\t$tcdb\t$name\t$keggcpd\t$chebcpd\t$hmdbcpd\t$pubccpd\t$inchcpd\t$bioccpd\n";
}
#######################################################	
#######################################################	






#######################################################	
#OUTPUT BIOCYC RXN DATA
#######################################################	
$time=localtime;
print "OUTPUT BioCyc rxn time $time\n";
print OUTBCR "rxn\tdb_src\trxn_dir\talt_rhea\talt_kegg\talt_bioc\t";
print OUTBCR "rhea_lcpd\trhea_lloc\trhea_ltrn\trhea_rcpd\trhea_rloc\trhea_rtrn\t";
print OUTBCR "kegg_lcpd\tkegg_lloc\tkegg_ltrn\tkegg_rcpd\tkegg_rloc\tkegg_rtrn\t";
print OUTBCR "bioc_lcpd\tbioc_lloc\tbioc_ltrn\tbioc_rcpd\tbioc_rloc\tbioc_rtrn\tUPIDs\tECs\n";
foreach my $rxn (sort(keys %RXN_ALTS)){
	$rhea_lcpd=''; 	$kegg_lcpd=''; 	$bioc_lcpd='';
	$rhea_lloc=''; 	$kegg_lloc=''; 	$bioc_lloc='';
	$rhea_ltrn=''; 	$kegg_ltrn=''; 	$bioc_ltrn='';
	$rhea_rcpd=''; 	$kegg_rcpd=''; 	$bioc_rcpd='';
	$rhea_rloc=''; 	$kegg_rloc=''; 	$bioc_rloc='';
	$rhea_rtrn=''; 	$kegg_rtrn=''; 	$bioc_rtrn='';
	$rrx='';	$krx='';	$brx='';
	$dir=$RXN_DIR{$rxn};
	$sorc=$RXN_ALTS{$rxn}{$rxn};


       foreach my $cpd (sort{$LRXN_CPD{$cpd}{$a}<=>$LRXN_CPD{$cpd}{$b} || $LRXN_CPD{$cpd}{$a} cmp $LRXN_CPD{$cpd}{$b}} keys %{$LRXN_CPD{$rxn}}){
                if($cpd !~/\w/){next;}
                if($TRANS_RXN{$rxn}{$cpd}=~/PORT/){ $tdir=$TRANS_RXN{$rxn}{$cpd};} else{$tdir="NOPORT";} #import/export/biport/noport
                if( $LRXN_CPD{$rxn}{$cpd}=~/SIDE/){ $tloc=$LRXN_CPD{$rxn}{$cpd};}  else{$tloc="INSIDE";} #inside or outside
                foreach my $alt (sort{$CMPD_ALTS{$cpd}{$a}<=>$CMPD_ALTS{$cpd}{$b} || $CMPD_ALTS{$cpd}{$a} cmp $CMPD_ALTS{$cpd}{$b}} keys %{$CMPD_ALTS{$cpd}}){
                           if($alt=~/^C\d+$/){          $kegg_lcpd.=$alt.";"; $kegg_lloc.=$tloc.";"; $kegg_ltrn.=$tdir.";";}
                        elsif($alt=~/^CHEBI.\d+$/){     $rhea_lcpd.=$alt.";"; $rhea_lloc.=$tloc.";"; $rhea_ltrn.=$tdir.";";}
                        elsif($alt=~/INCHI\:\w+|HMDB\d+|CID\:\d+/){next;}
                        else{                           $bioc_lcpd.=$alt.";"; $bioc_lloc.=$tloc.";"; $bioc_ltrn.=$tdir.";";}
                }
        }
        foreach my $cpd (sort{$RRXN_CPD{$cpd}{$a}<=>$RRXN_CPD{$cpd}{$b} || $RRXN_CPD{$cpd}{$a} cmp $RRXN_CPD{$cpd}{$b}} keys %{$RRXN_CPD{$rxn}}){
                if($cpd !~/\w/){next;}
                if($TRANS_RXN{$rxn}{$cpd}=~/PORT/){ $tdir=$TRANS_RXN{$rxn}{$cpd};}  else{$tdir="NOPORT";} #import/export/biport/noport
                if( $RRXN_CPD{$rxn}{$cpd} =~/SIDE/){ $tloc=$RRXN_CPD{$rxn}{$cpd}; } else{$tloc="INSIDE";} #inside or outside
                foreach my $alt (sort{$CMPD_ALTS{$cpd}{$a}<=>$CMPD_ALTS{$cpd}{$b} || $CMPD_ALTS{$cpd}{$a} cmp $CMPD_ALTS{$cpd}{$b}} keys %{$CMPD_ALTS{$cpd}}){
                           if($alt=~/^C\d+$/){          $kegg_rcpd.=$alt.";"; $kegg_rloc.=$tloc.";"; $kegg_rtrn.=$tdir.";";}
                        elsif($alt=~/^CHEBI.\d+$/){     $rhea_rcpd.=$alt.";"; $rhea_rloc.=$tloc.";"; $rhea_rtrn.=$tdir.";";}
                        elsif($alt=~/INCHI\:\w+|HMDB\d+|CID\:\d+/){next;}
                        else{                           $bioc_rcpd.=$alt.";"; $bioc_rloc.=$tloc.";"; $bioc_rtrn.=$tdir.";";}
                }
        }



	%ECS=(); %UPIDS=(); @UPIDS=(); $ec=''; $ex='';
	foreach my $ec   (keys %{$RXN_EC{$rxn}}){	$ECS{$ec}++; 
	  foreach my $upid (keys %{$EC2UPID{$ec}}){	$UPIDS{$upid}++;}}
	foreach my $upid (keys %{$RXN_UPID{$rxn}}){	$UPIDS{$upid}++;}
	foreach my $ex (sort{$ECS{$b}<=>$ECS{$a}} keys %ECS){ if($ex=~/\d/){$ec = $ex; last;} }
	foreach my $upid (keys %UPIDS){push(@UPIDS,$upid);}
	@UPIDS=nsort(@UPIDS);
	$upid=join(";",@UPIDS);
	foreach my $alt (sort(keys %{$RXN_ALTS{$rxn}})){
		   if($alt =~ /^\d+$/){ $rrx.=$alt.";"; }
		elsif($alt =~ /^R\d+$/){$krx.=$alt.";"; }
		elsif($alt =~ /RXN/){	$brx.=$alt.";"; }
		else{}
	}
	$one="$rxn\t$sorc\t$dir\t$rrx\t$krx\t$brx\t";
	$two="$rhea_lcpd\t$rhea_lloc\t$rhea_ltrn\t$rhea_rcpd\t$rhea_rloc\t$rhea_rtrn\t";
	$thr="$kegg_lcpd\t$kegg_lloc\t$kegg_ltrn\t$kegg_rcpd\t$kegg_rloc\t$kegg_rtrn\t";
	$for="$bioc_lcpd\t$bioc_lloc\t$bioc_ltrn\t$bioc_rcpd\t$bioc_rloc\t$bioc_rtrn\t$upid\t$ec";
	$out="$one$two$thr$for";
	$out=~s/\;+\t/\t/g;
	$out=~s/\;+$//g;
	print OUTBCR "$out\n";
}
undef(%CMPD_ALTS); undef(%CMPD_CHAR); undef(%CMPD_FORM); undef(%CMPD_MASS); undef(%CMPD_NAME); 
undef(%LRXN_CPD);  undef(%RRXN_CPD);  undef(%TRANS_CPD); undef(%LCPD_RXN);  undef(%RCPD_RXN);  undef(%TRANS_RXN); 
undef(%RXN_ALTS);  undef(%RXN_DIR);   undef(%RXN_EC); undef(%UPID_RXN); undef(%RXN_UPID);
###################################################################################
##############     DONE INPUT BIOCYC COMPOUNDS AND REACTIONS   ####################
###################################################################################










#####################################################################
########################   SUBROUTINES   ############################
#####################################################################
die;
#CHECK OUTPUTS CODE -- KEEP! - turn into subroutine?

#$cnt=0; $fkc  =keys %CMPD_FORM;	#foreach my $cpd (keys %CMPD_FORM){  						print "CMPD_FORM cpd $cpd val $CMPD_FORM{$cpd}\n";  $cnt++; if($cnt>100){last;}}
#$cnt=0; $mkc  =keys %CMPD_MASS;	#foreach my $cpd (keys %CMPD_MASS){  						print "CMPD_MASS cpd $cpd val $CMPD_MASS{$cpd}\n";  $cnt++; if($cnt>100){last;}}
#$cnt=0; $ckc  =keys %CMPD_CHAR;	#foreach my $cpd (keys %CMPD_CHAR){  						print "CMPD_CHAR cpd $cpd val $CMPD_CHAR{$cpd}\n";  $cnt++; if($cnt>100){last;}}
#$cnt=0; $akc  =keys %CMPD_ALTS;	#foreach my $cpd (keys %CMPD_ALTS){  foreach my $alt (keys %{$CMPD_ALTS{$cpd}}){print "CMPD_ALTS cpd $cpd alt $alt val $CMPD_ALTS{$cpd}{$alt}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $nkc  =keys %CMPD_NAME;	#foreach my $cpd (keys %CMPD_NAME){  foreach my $nam (keys %{$CMPD_NAME{$cpd}}){print "CMPD_NAME cpd $cpd nam $nam val $CMPD_NAME{$cpd}{$nam}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $lcx  =keys %LCPD_RXN;	#foreach my $cpd (keys %LCPD_RXN){   foreach my $rxn (keys %{$LCPD_RXN{$cpd}}){ print "LCPD_RXN  cpd $cpd rxn $rxn val $LCPD_RXN{$cpd}{$rxn}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $rcx  =keys %RCPD_RXN;	#foreach my $cpd (keys %RCPD_RXN){   foreach my $rxn (keys %{$RCPD_RXN{$cpd}}){ print "RCPD_RXN  cpd $cpd rxn $rxn val $RCPD_RXN{$cpd}{$rxn}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $trcpd=keys %TRANS_CPD; #foreach my $cpd (keys %TRANS_CPD){  foreach my $rxn (keys %{$TRANS_CPD{$cpd}}){print "TRANS_CPD cpd $cpd rxn $rxn val $TRANS_CPD{$cpd}{$rxn}\n"; } $cnt++; if($cnt>100){last;}}
##-
#$cnt=0; $rdir =keys %RXN_DIR;	#foreach my $rxn (keys %RXN_DIR){    						print "RXN_DIR   rxn $rxn val $RXN_DIR{$rxn}\n"; $cnt++; if($cnt>100){last;}}
#$cnt=0; $ralt =keys %RXN_ALTS;	#foreach my $rxn (keys %RXN_ALTS){   foreach my $alt (keys %{$RXN_ALTS{$rxn}}){ print "RXN_ALTS  rxn $rxn alt $alt val $RXN_ALTS{$rxn}{$alt}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $rcpd =keys %RRXN_CPD;	#foreach my $rxn (keys %RRXN_CPD){   foreach my $cpd (keys %{$RRXN_CPD{$rxn}}){ print "RRXN_CPD  rxn $rxn cpd $cpd val $RRXN_CPD{$rxn}{$cpd}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $lcpd =keys %LRXN_CPD;	#foreach my $rxn (keys %LRXN_CPD){   foreach my $cpd (keys %{$LRXN_CPD{$rxn}}){ print "LRXN_CPD  rxn $rxn cpd $cpd val $LRXN_CPD{$rxn}{$cpd}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $trrxn=keys %TRANS_RXN; #foreach my $rxn (keys %TRANS_RXN){  foreach my $cpd (keys %{$TRANS_RXN{$rxn}}){print "TRANS_RXN rxn $rxn cpd $cpd val $TRANS_RXN{$rxn}{$cpd}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $upid =keys %RXN_UPID;	#foreach my $rxn (keys %RXN_UPID){   foreach my $upx (keys %{$RXN_UPID{$rxn}}){ print "RXN_UPID rxn $rxn upx $upx val $RXN_UPID{$rxn}{$upx}\n"; } $cnt++; if($cnt>100){last;}}
#$cnt=0; $ec   =keys %RXN_EC;	#foreach my $rxn (keys %RXN_EC){     foreach my $ec (keys %{$RXN_EC{$rxn}}){ 	print "RXN_EC rxn $rxn ec $ec val $RXN_EC{$rxn}{$ec}\n"; } $cnt++; if($cnt>100){last;}}
#print "cpd info fkc $fkc mkc $mkc ckc $ckc akc $akc nkc $nkc lcx $lcx rcx $rcx trcpd $trcpd\n";
#print "rxn info rdir $rdir ralt $ralt lcpd $lcpd rcpd $rcpd trrxn $trrxn upid $upid ec $ec\n";

 


#FIX LOWQUAL CPD NAMES
sub CleanName{
	$nameX = $_[0];
        #remove junk punctuation/standardize
	$sta=0; $end=1;
	while($end ne $sta){
		$sta=$nameX;
		#swap greek symbols for text
		for my $g (0..$#GREEKL){ 	#fix pathbank and other greek symbols
			if($nameX =~/($GREEKS[$g])/){ 
				$nameX =~ s/$GREEKS[$g]/$GREEKL[$g]/g; 
		}	}
		$nameX =~ s/\%2B(\d*)/$1\+/g; 	#fix html +/- code (HMDB db)
		$nameX =~ s/\%2D(\d*)/$1\-/g; 	#fix html +/- code (HMDB db)
		
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
			#remove junk before name ex: 1_1_2_DI_5Z_8Z_11Z_14Z_EICOSATETRAENOYL_SN_GLYCERO_3_PHOSPHO_3_1_9Z_11Z_OCTADECADIENOYL_2_9Z_HEXADECENOYL_SN_GLYCERO_3_PHOSPHO_SN_GLYCEROL
		if($alt!~/[A-Z]{7,}/){$alt=$name; $alt=~s/[\W\_]+/\_/g;} #name too short after remove junk
	       	@PARTS = split("\_", $alt);
	       	$mp=0;
	       	foreach $p (@PARTS){ 
			$p=~/([A-Z]+)/; 
			if(length($1)>$mp){ $mp=length($1); } 
		}
		#$skip=0; foreach my $an (keys %LEN){if($an =~ /$alt/ && $an ne $alt){ $skip=1; last;}}
		#if($skip==1){next;}
	       	$LEN{$alt}=$mp;	#whats longest letter stretch
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
		elsif(length($name) >= 50){		    $LONG{$name} =1; }
		else{ print "badname $name\n"; }
	}


	#COMBINE NESTED NAMES
	@NAMES=();
	foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %GOOD){
		$do=1;
		for my $i (0..$#NAMES){ #loop through added names, remove nested
			if($NAMES[$i]=~/$name/){$do=0; last;} 			#the name is contained in another name, skip/redundant
			if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
		}
		if($do==1){push(@NAMES,$name);} 				#no name containment
 		if($NAMES[9]=~/\w/){last;}					#limit to 10 good names
	}
	foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %LONG){
		$do=1;
		for my $i (0..$#NAMES){						#excluding nested names
			if($NAMES[$i]=~/$name/){$do=0; last;} 			#the name is contained in another name, skip/redundant
			if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
		}
		if($do==1){push(@NAMES,$name);} 				#no name containment 
		if($NAMES[10]=~/\w/){last;} 					#if 10 "good" names use only 1 long
	}
	foreach my $name (sort{ $ODD{$a}<=>$ODD{$b} || $LEN{$b}<=>$LEN{$a} } keys %SHORT){
		$do=1;
		for my $i (0..$#NAMES){						#excluding nested names
			if($NAMES[$i]=~/$name/){$do=0; last;} 			#the name is contained in another name, skip/redundant
			if($name=~/$NAMES[$i]/){$do=0; $NAMES[$i]=$name; last;} #the name contains a prior name, replace with name-longer
		}
		if($do==1){push(@NAMES,$name);} 				#no name containment 
		if($NAMES[11]=~/\w/){last;}                                      #if 10 "good" and 1 long names use only 1 short
	}
	@NAMES=nsort(@NAMES);
       	return(@NAMES);
}
#another stupid name: <synonym>alpha,alpha,Alpha',alpha'-tetramethyl-5-(1H-1,2,4-triazol-1-ylmethyl)-m-benzenediacetonitrile</synonym>
#or -3-5E_7Z_11Z_14Z-9-10-BENZOYLOXY-1_2_6A_6B_9_9_12A-HEPTAMETHYL-13-OXO-1_2_3_4_4A_5_6_6A_6B_7_8_8A_9_10_11_12_12A_12B_13_14B-ICOSAHYDROPICENE-4A-CARBOXYLATE
#or people cant spell:
#good DIMETHYLPHOSPHOROAMIDOTHIOATE
#good DIMETHYL_PHOSPHORAMIDOTHIOLATE
#good DIMETHYL_PHOSPHOROAMIDOTHIOATE
#good DIMETHYL_PHOSPHORAMIDOTHIOATE
#good PHOSPHORAMIDOTHIOATE_O_S_DIMETHYL_ESTER
#good METHYL_PHOSPHORAMIDOTHIOATE
#good DIMETHYL_PHOSPHORAMIDOTHIOLIC_ACID
#good METHYL_PHOSPHORAMIDOTHIOATE_MEO_MES
#good PHOSPHORAMIDOTHIOIC_ACID_O_S_DIMETHYL_ESTER
#good METHYL_PHOSPHORAMIDOTHIOIC_ACID
#good DIMETHYL_AMIDOTHIOPHOSPHORIC_ACID
#good THIOPHOSPHORAMIDATE_O_S_DIMETHYL_ESTER
#good DIMETHYL_AMIDOTHIOPHOSPHATE
#good THIOPHOSPHORAMIDIC_ACID_O_S_DIMETHYL_ESTER
#OR JUST STUPID LONG:
#CARBOXY_6_CARBOXY_DIHYDROXY_HYDROXY_HYDROXYOCTADEC_9_EN_1_YLIDENE_AMINO_OCTADEC_4_EN_1_YL_OXY_2_HYDROXYMETHYL_OXAN_3_YL_OXY_3_HYDROXY_HYDROXY_6_HYDROXYMETHYL_OXOPROPYL_TRIHYDROXY_6_HYDROXYMETHYL_OXAN_2_YL_OXY_OXAN_2_YL_OXY_6_HYDROXYMETHYL_OXAN_4_YL_OXY_4_HYDROXY_HYDROXYETHYLIDENE_AMINO_OXAN_2_YL_DIHYDROXYPROPAN_2_YL_OXY_4_HYDROXY_HYDROXYETHYLIDENE_AMINO_OXAN_2_YL_DIHYDROXYPROPAN_2_YL_OXY_4_HYDROXY_HYDROXYETHYLIDENE_AMINO_TRIHYDROXYPROPYL_OXANE_2_CARBOXYLATE

#####################################################################
########################   END SUBROUTINES   ########################
#####################################################################





__END__











__END__
#TO get remaining pubchem info:
#1. go to URLs below, select "send to" then "file" then "properties (text)" then "create file"
#2. pause or cancel download (in chrome CTRL+J) then right click and copy the LONG url
#3. download in a shell: wget -O [whichever database name] "https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=pccompound...[the LONG URL]..."
#yes, this is annoying and manual, you can try downloading tons of HUGE xml files from pubchem's ftp, but that will take even longer!
#chembl and chembridge both have 2M+ compounds - so you may not need the others - but those two have distinct compounds, and the others might as well
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22CHEBI%22%5BSourceName%5D
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22KEGG%22%5BSourceName%5D
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22BIOCYC%22%5BSourceName%5D
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22Human+Metabolome+Database+(HMDB)%22%5BSourceName%5D
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22ChEMBL%22%5BSourceName%5D
#https://www.ncbi.nlm.nih.gov/pccompound/?term=%22ChemBridge%22%5BSourceName%5D
#stop
#the below wgets are an example, they will not work for you - sorry you have to manual download, change the -O file_names.txt
#qx{wget -O pubchem_biocyc.txt 	"https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=pccompound&HistoryId=MCID_6245bfb19b12c849df7029a0&QueryKey=6&Sort=CIDA&Filter=all&CompleteResultCount=25449&Mode=file&View=property&p$l=Email&portalSnapshot=%2Fprojects%2FPubChem%2FPubChem_dbs%401.20&BaseUrl=&PortName=live&FileName="};
#qx{wget -O pubchem_cbridge.txt "https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=pccompound&HistoryId=MCID_6245cb6ef0bea96738289a33&QueryKey=1&Sort=&Filter=all&CompleteResultCount=1734455&Mode=file&View=property&p$l=Email&portalSnapshot=%2Fprojects%2FPubChem%2FPubChem_dbs%401.20&BaseUrl=&PortName=live&FileName="};
#qx{wget -O pubchem_chembl.txt 	"https://www.ncbi.nlm.nih.gov/portal/utils/file_backend.cgi?Db=pccompound&HistoryId=MCID_6245cb6ef0bea96738289a33&QueryKey=1&Sort=CIDA&Filter=all&CompleteResultCount=2136922&Mode=file&View=property&p$l=Email&portalSnapshot=%2Fprojects%2FPubChem%2FPubChem_dbs%401.20&BaseUrl=&PortName=live&FileName="};



