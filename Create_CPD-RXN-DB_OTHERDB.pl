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
qx{wget -N https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip};    qx{unzip -u pathbank_all_metabolites.csv.zip};
qx{wget -N https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip};      qx{unzip -u hmdb_metabolites.zip};
qx{wget -N https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl};
qx{wget -O CID_cheb_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22chebi%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_chebi%22}'};
qx{wget -O CID_bioc_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22biocyc%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_biocyc%22}'};
qx{wget -O CID_hmdb_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22hmdb%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_hmdb%22}'};
qx{wget -O CID_kegg_summary.csv 'https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22substance%22,%22where%22:{%22ands%22:[{%22*%22:%22kegg%22}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_substance_text_kegg%22}'};
qx{wget -N https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz};	   qx{gunzip -fk CID-InChI-Key.gz};
qx{wget -N https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Mass.gz};		   qx{gunzip -fk CID-Mass.gz};
qx{wget -N https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz};  qx{gunzip -fk CID-Synonym-filtered.gz};


#OPEN/CHECK INPUTS:
$time=localtime;
print "OPEN / CHECK INPUTS time $time\n";
#universal
open(INTRCH,	"getSubstrates.py")		||die "unable to open getSubstrates.py: $!\n";
$ec_up = `wget -q -O - https://ftp.expasy.org/databases/enzyme/enzyme.dat`;
#specific
open(INPBANK, 	"pathbank_all_metabolites.csv")	||die "unable to open pathbank_all_metabolites.csv: $!\n";
open(INHMDB, 	"hmdb_metabolites.xml")		||die "unable to open hmdb_metabolites.xml: $!\n";
open(INPUBBIOC, "CID-biocyc_summary.csv")	||die "Unable to open CID_biocyc_summary.csv\n";
open(INPUBCHEB, "CID-chebi_summary.csv") 	||die "Unable to open CID_chebi_summary.csv\n";
open(INPUBKEGG, "CID-kegg_summary.csv")	 	||die "Unable to open CID_kegg_summary.csv\n";
open(INPUBHMDB, "CID-hmdb_summary.csv")	 	||die "Unable to open CID_hmdb_summary.csv\n";
open(INPUBINCH, "CID-InChI-Key")	 	||die "Unable to open CID-InChI-Key\n";
open(INPUBMASS, "CID-Mass")		 	||die "Unable to open CID-Mass\n";
open(INPUBSYN, 	"CID-Synonym-filtered")	 	||die "Unable to open CID-Synonym-filtered\n";

#OPEN/CHECK OUTPUTS:
open(OUTPBC, ">", "PATHBANK_CPD_DB.txt")||die;
open(OUTHMC, ">", "HMDB_CPD_DB.txt")||die;
open(OUTPCC, ">", "PUBCHEM_CPD_DB.txt")||die;
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
#distribute upids to oboletes
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
######################   INPUT PATHBANK COMPOUNDS	 #########################
##################################################################################
# DO PATHBANK FIRST - least reliable, can be replaced with other DBs
# ALL OTHERs compounds: KEGG,CHEBI,PUBCHEM,HMDB,BIOCYC - will be replaced with database-specific information
# PLUS pathbank such a pain to parse, less trusty

$on=0;
$time=localtime;
print "INPUT pathbank_all_metabolites.csv time $time\n";
while(<INPBANK>){
	if($_!~/\w/){next;}
	if($_=~/^PATHBANK.ID/i){next;}
	if($on%200000==0){$time=localtime; print "on $on time $time\n";} $on++;
 	$_=~s/[\r\n]+//;
	$_=uc($_);
	$line=$_;
	$line=~s/.*\,PW_C\d+\,//;

	#swap greek symbols for text
	for my $g (0..$#GL){ 
		if($line =~/($GS[$g])/){ 
			$line =~ s/$GS[$g]/$GL[$g]/g; 
	}	}

	#clear out empty "" from csv mess with parse
	while($line=~/\,\"\"\,/){$line=~s/\,\"\"\,/\,\,/;}

	#get formula before changing cell formatting
	$form1=''; if($line =~ /INCHI.[^\/]+\/([^\/]+)\//){ $form1=$1; $form1=~s/\W+/\_/g; $form1=CleanNames($form1); }

	#fix commas and non-word from cell values
	while($line=~/(\"[^\"]+\")/){ 
		$found=$1;
		$new=$found; 
		$new=~s/\"//g;
		$new=~s/[^\w+\-\+]+/\_/g;
		$line=~s/\Q$found\E/$new/;
	}
	while($line=~/[^\,\w\-\+]+/){$line=~s/[^\,\w\-\+]+/\_/;}
	$line=~s/\_+/\_/g;

	%ALTS=();
	@LINES=();	
	@LINES = split(",", $line);
	if($#LINES !=10 || $LINES[0]=~/ELECTRON/){next;} #electron containing - all have wrong collumn numbers
	$name=''; if($LINES[0]=~/\w/){		  $name=$LINES[0]; $name=CleanNames($name);}
	$form=''; if($LINES[6]=~/\w/){
			$form2=$LINES[6]; 
			$form2=CleanNames($form2);
			if($form !~ /\w/){$form=$form2;}
			if($form1 =~/\w/ && $form !~ /\w/){$form=$form1;}
		  }

	$id='';   if($LINES[1]=~/(HMDB\d+)/){	  $id=$1; while(length($id)<11){$id=~s/HMDB/HMDB0/;} 	$ALTS{$id}=1;}
	$id='';   if($LINES[2]=~/(C\d\d\d\d\d)/){ $id=$1;						$ALTS{$id}=1;}
	$id='';   if($LINES[3]=~/^(\d+)$/){	  $id=$1; $id="CHEBI:".$id; 				$ALTS{$id}=1;}
	$id='';	  if($LINES[10]=~/([A-Z\-]+)/){	  $id=$1; $id="INCHI:".$id; 				$ALTS{$id}=1;}

	#pathbank is general db so no alts - all are primary ids
	foreach my $id (keys %ALTS){
		$id =~ s/\s+//g; 
		if($id !~ /\w/){next;}
		$CMPD_ALTS{$id}{$id}="PB";
		if($name=~/\w/){$CMPD_NAME{$id}{$name}="PB";}
		   if($form =~/\w/){$CMPD_FORM{$id}=$form;}
		elsif($form1=~/\w/){$CMPD_FORM{$id}=$form1;}
		else{$form='';}
		foreach my $id2 (keys %ALTS){
			$id2 =~ s/\s+//g;
			if($id2 !~ /\w/){next;}
			$CMPD_ALTS{$id}{$id2}="PB";
}	}	}	


#OUTPUT PATHBANK - DO NOT MIX THIS OUTPUT SCRIPT PIECE WITH THE OTHERS
$time=localtime;
print "OUTPUT pathbank_all_metabolites.csv time $time\n";
print OUTPBC "cpd\tformula\tmass\tcharge\tdb_src\ttcdbs\tnames\tkeggcpd\tchebcpd\thmdbcpd\tpubccpd\tinchcpd\tbioccpd\n";
foreach my $cpd (sort(keys %CMPD_ALTS)){
	$fe=''; $fa=''; $form=''; $char=''; $mass=''; $name=''; $dupl=0;
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

	#sort identifiers by database
	foreach my $alt (sort(keys %{$CMPD_ALTS{$cpd}})){
		   if($alt=~/^C\d+$/){		if($keggcpd!~/$alt/){$keggcpd.=$alt.";"; }}
		elsif($alt=~/^CHEBI.\d+$/){	if($chebcpd!~/$alt/){$chebcpd.=$alt.";"; }}
		elsif($alt=~/^HMDB\d+$/){	if($hmdbcpd!~/$alt/){$hmdbcpd.=$alt.";"; }}
		elsif($alt=~/^CID.\d+$/){	if($pubccpd!~/$alt/){$pubccpd.=$alt.";"; }}
		elsif($alt=~/^INCHI.\S+$/){	if($inchcpd!~/$alt/){$inchcpd.=$alt.";"; }}
		else{				if($bioccpd!~/$alt/){$bioccpd.=$alt.";"; }}
	}
	#check for multiple identfiers
	if($keggcpd=~/\;.*\;/ || $chebcpd=~/\;.*\;/ || $hmdbcpd=~/\;.*\;/ || $pubccpd=~/\;.*\;/ || $inchcpd=~/\;.*\;/ || $bioccpd=~/\;.*\;/){$dupl=1;}
	if($dupl==1 && $form=~/\w/){
		#THERE IS NO PRIMARY CPDID FOR PATHBANK
		#USING FORMULA COUNT TO CHOOSE MOST LIKELY
		#GET EVERY ALT CPD FORMULA ELEMENTS
		%FORMS=(); %FOCO=(); %FOMO=(); $af=''; $bfo='';
		foreach my $alt (keys %{$CMPD_ALTS{$cpd}}){
			$alf=$CMPD_FORM{$alt};
			if($alf !~ /\w/){next;}
			$FORMS{$alf}=$alf;
			$af=$alf;
			$af=~s/[^A-Z]+//g;
			$FOCO{$af}++;
			$FOMO{$af}=$alt;
		}

		#sort decending max count for each formula keep top
		foreach my $af (sort{$FOCO{$b}<=>$FOCO{$a}} keys %FOCO){
			if($bfo eq ''){ $bfo=$af; }
			if($af ne $bfo){
				$alt=$FOMO{$af};
				$frm=$CMPD_FORM{$alt};
				delete($CMPD_ALTS{$cpd}{$alt});
				delete($FORMS{$frm});
				   if($alt=~/^C\d+$/){		$keggcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				elsif($alt=~/^CHEBI.\d+$/){	$chebcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				elsif($alt=~/^HMDB\d+$/){	$hmdbcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				elsif($alt=~/^CID.\d+$/){	$pubccpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				elsif($alt=~/^INCHI.\S+$/){	$inchcpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				else{				$bioccpd=~s/^$alt\;|\;$alt\;|\;$alt$|^$alt$//;}
				next;
		}	}

		#if same elements but different number replace with "R"
		$form='';
		@ATOMS=split("",$bfo);
		foreach my $atom (@ATOMS){
			%CNT=();
			foreach my $f (keys %FORMS){
				   if($FORMS{$f} =~ s/^$atom(\d+)//){	$CNT{$1}=1;}
				elsif($FORMS{$f} =~ s/^$atom//){ 	$CNT{""}=1;}
				else{ delete($FORMS{$f});}
			}
			if(keys %CNT == 1){   foreach my $c (keys %CNT){$form.=$atom.$c;}}
			else{ 						$form.=$atom."R";}
	}	}								
	$keggcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $keggcpd); @CPDS=nsort(@CPDS); $keggcpd=join(";",@CPDS);
	$chebcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $chebcpd); @CPDS=nsort(@CPDS); $chebcpd=join(";",@CPDS);
	$cc=@CPDS; %TCDBS=(); @TCDBS=();  $tcdb='';
	if($cc==1){
		foreach my $ch (@CPDS){ foreach my $tc (sort(keys %{$CHEB_TCDB{$ch}})){ $TCDBS{$tc}=1; }}
		if(keys %TCDBS > 0){ foreach my $tc (sort(keys %TCDBS)){ push(@TCDBS,$tc); } @TCDBS=nsort(@TCDBS); $tcdb=join(";",@TCDBS);}
	}
	$hmdbcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $hmdbcpd); @CPDS=nsort(@CPDS); $hmdbcpd=join(";",@CPDS);
        $pubccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $pubccpd); @CPDS=nsort(@CPDS); $pubccpd=join(";",@CPDS);
	$inchcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $inchcpd); @CPDS=nsort(@CPDS); $inchcpd=join(";",@CPDS);
	$bioccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $bioccpd); @CPDS=nsort(@CPDS); $bioccpd=join(";",@CPDS);
	print OUTPBC "$cpd\t$form\t$mass\t$char\t$sorc\t$tcdb\t$name\t$keggcpd\t$chebcpd\t$hmdbcpd\t$pubccpd\t$inchcpd\t$bioccpd\n";
}

### ISSUES WITH PATHBANK:
#pathbank horrible parsing mess
#CANT TRACK ++ for names/alts/cpds because pathbank has huge duplicate issue
#Using database abbreviation instead
#another error in pathbank - shows how same kegg and chebi cpd IDs are associated with different formulas
#in this and many cases, the number of different formulas/mixed ids can be substantial
#SMP0027679,Metabolic,Homo sapiens,PW_C063167,HMDB0116024	,C00269,17962,,,C52H85N3O15P2,"{[(2R,3R,5R)-5-(4-amino-2-oxo-1,2-dihydropyrimidin-1-yl)-3,4-dihydroxyoxolan-2-yl]methoxy}({[(2R)-3-[(7Z,10Z,13Z,16Z,19Z)-docosa-7,10,13,16,19-pentaenoyloxy]-2-[(9Z)-octadec-9-enoyloxy]propoxy](hydroxy)phosphoryl}oxy)phosphinic acid",[H][C@@](COC(=O)CCCCC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CC)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C(O)[C@H]1O)N1C=CC(N)=NC1=O)OC(=O)CCCCCCC\C=C/CCCCCCCC,"InChI=1S/C52H85N3O15P2/c1-3-5-7-9-11-13-15-17-19-20-21-22-24-25-27-29-31-33-35-37-47(56)65-41-44(68-48(57)38-36-34-32-30-28-26-23-18-16-14-12-10-8-6-4-2)42-66-71(61,62)70-72(63,64)67-43-45-49(58)50(59)51(69-45)55-40-39-46(53)54-52(55)60/h5,7,11,13,17-19,21-23,25,27,39-40,44-45,49-51,58-59H,3-4,6,8-10,12,14-16,20,24,26,28-38,41-43H2,1-2H3,(H,61,62)(H,63,64)(H2,53,54,60)/b7-5-,13-11-,19-17-,22-21-,23-18-,27-25-/t44-,45-,49+,50?,51-/m1/s1",SKIDWIJLVZVLQC-SCINYDAUSA-N
#SMP0027680,Metabolic,Homo sapiens,PW_C0455",			,C00269,17962,,,C46H81N3O15P2,"{[(2R,3R,5R)-5-(4-amino-2-oxo-1,2-dihydropyrimidin-1-yl)-3,4-dihydroxyoxolan-2-yl]methoxy}({[(2R)-3-[(9Z)-hexadec-9-enoyloxy]-2-[(11Z)-octadec-11-enoyloxy]propoxy](hydroxy)phosphoryl}oxy)phosphinic acid",[H][C@@](COC(=O)CCCCCCC\C=C/CCCCCC)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C(O)[C@H]1O)N1C=CC(N)=NC1=O)OC(=O)CCCCCCCCC\C=C/CCCCCC,"InChI=1S/C46H81N3O15P2/c1-3-5-7-9-11-13-15-17-18-20-22-24-26-28-30-32-42(51)62-38(35-59-41(50)31-29-27-25-23-21-19-16-14-12-10-8-6-4-2)36-60-65(55,56)64-66(57,58)61-37-39-43(52)44(53)45(63-39)49-34-33-40(47)48-46(49)54/h13-16,33-34,38-39,43-45,52-53H,3-12,17-32,35-37H2,1-2H3,(H,55,56)(H,57,58)(H2,47,48,54)/b15-13-,16-14-/t38-,39-,43+,44?,45-/m1/s1",YOVNZDXBUXIJEK-HOXHGPEGSA-N
#SMP0027680,Metabolic,Homo sapiens,PW_C063170,HMDB0116027	,C00269,17962,,,C56H85N3O15P2,"{[(2R,3R,5R)-5-(4-amino-2-oxo-1,2-dihydropyrimidin-1-yl)-3,4-dihydroxyoxolan-2-yl]methoxy}({[(2R)-2-[(4Z,7Z,10Z,13Z,16Z)-docosa-4,7,10,13,16-pentaenoyloxy]-3-[(7Z,10Z,13Z,16Z,19Z)-docosa-7,10,13,16,19-pentaenoyloxy]propoxy](hydroxy)phosphoryl}oxy)phosphinic acid",[H][C@@](COC(=O)CCCCC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CC)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C(O)[C@H]1O)N1C=CC(N)=NC1=O)OC(=O)CC\C=C/C\C=C/C\C=C/C\C=C/C\C=C/CCCCC,"InChI=1S/C56H85N3O15P2/c1-3-5-7-9-11-13-15-17-19-21-23-25-27-29-31-33-35-37-39-41-51(60)69-45-48(72-52(61)42-40-38-36-34-32-30-28-26-24-22-20-18-16-14-12-10-8-6-4-2)46-70-75(65,66)74-76(67,68)71-47-49-53(62)54(63)55(73-49)59-44-43-50(57)58-56(59)64/h5,7,11-14,17-20,23-26,29-32,36,38,43-44,48-49,53-55,62-63H,3-4,6,8-10,15-16,21-22,27-28,33-35,37,39-42,45-47H2,1-2H3,(H,65,66)(H,67,68)(H2,57,58,64)/b7-5-,13-11-,14-12-,19-17-,20-18-,25-23-,26-24-,31-29-,32-30-,38-36-/t48-,49-,53+,54?,55-/m1/s1",ZLQYPSZULACELO-GZNMQYDESA-N
undef(%CMPD_ALTS); undef(%CMPD_CHAR); undef(%CMPD_FORM); undef(%CMPD_MASS); undef(%CMPD_NAME); 
undef(%LRXN_CPD);  undef(%RRXN_CPD);  undef(%TRANS_CPD); undef(%LCPD_RXN);  undef(%RCPD_RXN);  undef(%TRANS_RXN); 
undef(%RXN_ALTS);  undef(%RXN_DIR);   undef(%RXN_EC); undef(%UPID_RXN); undef(%RXN_UPID);
##################################################################################
######################   DONE INPUT PATHBANK COMPOUNDS	 #########################
##################################################################################











#################################################################################
####################   INPUT HMDB COMPOUNDS AND REACTIONS   #####################
#################################################################################
$on=0;
$time=localtime;
print "INPUT hmdb_metabolites.xml time $time\n";
while(<INHMDB>){
        if($_!~/\w/){next;}
        $_=~s/[\r\n]+//;
        $_=uc($_);
	if($on%10000000==0){$time=localtime; print "on $on time $time\n";} $on++;
	if($_=~/\<METABOLITE\>/i){		$inmet=1; %ALTS=(); %NA=(); $id=''; $alt=''; $mass=''; $name=''; $form=''; $char=''; $inchrg=0; $insec=0; 	next;}
	if($_=~/\<secondary.accessions/i){ 	$insec=1; }
	if($_=~/\<\/secondary.accessions/i){ 	$insec=0; }
	if($insec==1){ if($_=~/\<accession\>(HMDB\d+)/i){          	$alt=$1; while(length($alt)<11){$alt=~s/HMDB/HMDB0/;}  	$ALTS{$alt}=1;	next;}}
	if($inmet==1){
			    if($_=~/\<accession\>(HMDB\d+)/i){ 		$id=$1;  while(length($id)<11){$id=~s/HMDB/HMDB0/;} 	$ALTS{$id}=1;   next;}
				    if($_=~/\<inchikey\>([A-Z\-]+)/i){ 	$alt="INCHI:".$1;					$ALTS{$alt}=1; 	next;}
				       if($_=~/\<kegg.id\>(C\d{5})/i){ 	$alt=$1;						$ALTS{$alt}=1; 	next;}
			      if($_=~/\<pubchem.compound.id\>(\d+)/i){ 	$alt="CID:".$1; 					$ALTS{$alt}=1; 	next;}
				         if($_=~/\<chebi.id\>(\d+)/i){ 	$alt="CHEBI:".$1;					$ALTS{$alt}=1; 	next;}
	                              if($_=~/\<biocyc.id\>(.+?)\</i){ 	$alt=$1; 	  					$ALTS{$alt}=1; 	next;}
				    if($_=~/^\s{2}\<name\>([^\<]+)/i){ 	$name=$1; $name=CleanNames($name);			$NA{$name}=1;	next;}
				 if($_=~/^\s{4}\<synonym\>([^\<]+)/i){ 	$name=$1; $name=CleanNames($name);			$NA{$name}=1;   next;}
			       if($_=~/\<chemical.formula\>(.+?)\</i){ 	$form=$1; $form=CleanNames($form); 					next;}
	             if($_=~/\<average.molecular.weight\>([\d\.]+)/i){ 	$mass=$1;								next;}
		            if($_=~/\<kind\>physiological.charge\</i){ 	$inchrg=1; 								next;}
		      if($inchrg==1 && $_=~/\<value\>([\+\-]*)(\d+)/i){ $inchrg=0; $s=''; $s=$1;$n=''; $n=$2; $char=''; 
									if($n!=0 && $s ne "-"){$char=$n."+";} else{$char=$n.$s;}		next;}
	}
	if($_=~/\<\/METABOLITE\>/i){
		$inmet=0;
		$id =~ s/\s+//g; 
		if($id !~ /\w/){next;}
		$CMPD_ALTS{$id}{$id}="HM";
		foreach my $alt (keys %ALTS){ 
			$alt =~ s/\s+//g; if($alt !~ /\w/){next;}
			$CMPD_ALTS{$id}{$alt}="HM"; 
		}
		foreach my $nm (keys %NA){
		if($nm=~/\w/){	$CMPD_NAME{$id}{$nm}="HM";}}
		if($mass=~/\w/){$CMPD_MASS{$id}=$mass;}
		if($char=~/\w/){$CMPD_CHAR{$id}=$char;}
		if($form=~/\w/){$CMPD_FORM{$id}=$form;}
		next;
	}	
}


#OUTPUT
$time=localtime;
print "OUTPUT hmdb_metabolites.xml time $time\n";
print OUTHMC "cpd\tformula\tmass\tcharge\tdb_src\ttcdbs\tnames\tkeggcpd\tchebcpd\thmdbcpd\tpubccpd\tinchcpd\tbioccpd\n";
foreach my $cpd (sort(keys %CMPD_ALTS)){
	$fe=''; $fa=''; $form=''; $char=''; $mass=''; $name='';
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
	if(keys %TCDBS > 0){ foreach my $tc (sort(keys %TCDBS)){ push(@TCDBS,$tc); } @TCDBS=nsort(@TCDBS); $tcdb=join(";",@TCDBS);}}
	$hmdbcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $hmdbcpd); @CPDS=nsort(@CPDS); $hmdbcpd=join(";",@CPDS);
        $pubccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $pubccpd); @CPDS=nsort(@CPDS); $pubccpd=join(";",@CPDS);
	$inchcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $inchcpd); @CPDS=nsort(@CPDS); $inchcpd=join(";",@CPDS);
	$bioccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $bioccpd); @CPDS=nsort(@CPDS); $bioccpd=join(";",@CPDS);
	print OUTHMC "$cpd\t$form\t$mass\t$char\t$sorc\t$tcdb\t$name\t$keggcpd\t$chebcpd\t$hmdbcpd\t$pubccpd\t$inchcpd\t$bioccpd\n";
}
undef(%CMPD_ALTS); undef(%CMPD_CHAR); undef(%CMPD_FORM); undef(%CMPD_MASS); undef(%CMPD_NAME); 
undef(%LRXN_CPD);  undef(%RRXN_CPD);  undef(%TRANS_CPD); undef(%LCPD_RXN);  undef(%RCPD_RXN);  undef(%TRANS_RXN); 
undef(%RXN_ALTS);  undef(%RXN_DIR);   undef(%RXN_EC); undef(%UPID_RXN); undef(%RXN_UPID);
#################################################################################
####################   DONE HMDB COMPOUNDS AND REACTIONS   ######################
#################################################################################
 










#################################################################################
#####################	 INPUT PUBCHEM COMPOUNDS	#########################
#################################################################################
$time=localtime;
print "GET PUBCHEM COMPOUNDS time $time\n";
print "INPUT PubChem_substance_text_biocyc_summary.csv\n";
while(<INPUBBIOC>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	$_=~s/\"\"//g;
	$_ =~ /^[\s\d]*\,(\d+)\,BioCyc\,([^\,]+)\,/i;
	$cid="CID:".$1;
	$bid=$2;
        $cid="CID:".$1;
        $bid=$2;
	if($cid=~/\d/ && $bid=~/\d/){$CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; }
}

print "INPUT PubChem_substance_text_chebi_summary.csv\n";
while(<INPUBCHEB>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	$_=~s/\"\"//g;
	$_ =~ /^[\s\d]*\,(\d+)\,.*?(CHEBI\:\d+)/i;
        $cid="CID:".$1;
        $bid=$2;
	if($cid=~/\d/ && $bid=~/\d/){$CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; }
}

print "INPUT PubChem_substance_text_kegg_summary.csv\n";
while(<INPUBKEGG>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	$_=~s/\"\"//g;
	$_ =~ /^[\s\d]*\,(\d+)\,KEGG.*?(C\d\d\d\d\d)/i;
	$cid="CID:".$1;
	$bid=$2;
	if($cid=~/\d/ && $bid=~/\d/){$CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; }
}

print "INPUT PubChem_substance_text_hmdb_summary.csv\n";
while(<INPUBHMDB>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	$_=~s/\"\"//g;
	$_ =~ /^[\s\d]*\,(\d+)\,.*?HMDB.*?(HMDB\d+)/i;
	$cid="CID:".$1;
	$bid=$2;
	if($bid !~/HMDB/){next;} 
	while(length($bid)<11){$bid=~s/HMDB/HMDB0/;}
	if($cid=~/\d/ && $bid=~/\d/){$CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; }
}

print "INPUT CID-InChI-Key.gz\n";
while(<INPUBINCH>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	(my $cid, my $junk, my $bid)=split("\t",$_);
	if($cid!~/^\d+$/){next;}
	$cid="CID:".$cid;
	$bid="INCHI:".$bid;
	if($bid !~/[A-Z\-]/){next;} 
	if($cid=~/\d/ && $bid=~/\w/){$CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; }
}

print "INPUT CID-Mass.gz\n";
while(<INPUBMASS>){
	if($_!~/\w/){next;} 
	$_=~s/[\r\n]+//; 
	$_=uc($_); 
	(my $cid, my $form, my $mass, my $junk)=split("\t",$_); 
	if($cid!~/^\d+$/){next;}
	$cid="CID:".$cid;
	if($form=~/\w/){$CMPD_FORM{$cid}=CleanNames($form);;} 
	if($mass=~/\d/){$CMPD_MASS{$cid}=$mass;}
}

print "INPUT CID-Synonym-filtered.gz\n";
while(<INPUBSYN>){
	if($_!~/\w/){next;}
	$_=~s/[\r\n]+//;
	$_=uc($_);
	(my $cid, my $bid)=split("\t",$_);
	if($cid!~/^\d+$/){next;}
	$cid="CID:".$cid;
	if($bid=~/^CHEBI\:\d+$/){ $CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; next;}
	if($bid=~/^C\d{5}$/){ 	  $CMPD_ALTS{$cid}{$bid}="PC"; $CMPD_ALTS{$cid}{$cid}="PC"; next;}
	$xid="INCHI:".$bid; if(exists($CMPD_ALTS{$cid}{$xid})){next;} #get rid of inchikey from names
	if($bid=~/[A-Z]{8}/){ $name=CleanNames($bid); $CMPD_NAME{$cid}{$name}="PC";}

}
foreach my $cid (sort(keys %CMPD_NAME)){
	@NAMES=();
	foreach my $name (keys %{$CMPD_NAME{$cid}}){ push(@NAMES,$name);}
	@NAMES=BestName(@NAMES);
	delete($CMPD_NAME{$cid});
	foreach my $name (@NAMES){ $CMPD_NAME{$cid}{$name}="PC"; }
}


#OUTPUT PUBCHEM CPDS
$time=localtime;
print "OUTPUT PubChem time $time\n";
print OUTPCC "cpd\tformula\tmass\tcharge\tdb_src\ttcdbs\tnames\tkeggcpd\tchebcpd\thmdbcpd\tpubccpd\tinchcpd\tbioccpd\n";
foreach my $cpd (sort(keys %CMPD_ALTS)){
	$fe=''; $fa=''; $form=''; $char=''; $mass=''; $name='';
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
			#print "cpd $cpd form $form alf $alf alt $alt not match =delete\n";
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
	      foreach my $tc (sort(keys %TCDBS)){ push(@TCDBS,$tc); } @TCDBS=nsort(@TCDBS); $tcdb=join(";",@TCDBS);}
	}
	$hmdbcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $hmdbcpd); @CPDS=nsort(@CPDS); $hmdbcpd=join(";",@CPDS);
        $pubccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $pubccpd); @CPDS=nsort(@CPDS); $pubccpd=join(";",@CPDS);
	$inchcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $inchcpd); @CPDS=nsort(@CPDS); $inchcpd=join(";",@CPDS);
	$bioccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $bioccpd); @CPDS=nsort(@CPDS); $bioccpd=join(";",@CPDS);
	print OUTPCC "$cpd\t$form\t$mass\t$char\t$sorc\t$tcdb\t$name\t$keggcpd\t$chebcpd\t$hmdbcpd\t$pubccpd\t$inchcpd\t$bioccpd\n";
}
undef(%CMPD_ALTS); undef(%CMPD_CHAR); undef(%CMPD_FORM); undef(%CMPD_MASS); undef(%CMPD_NAME); 
undef(%LRXN_CPD);  undef(%RRXN_CPD);  undef(%TRANS_CPD); undef(%LCPD_RXN);  undef(%RCPD_RXN);  undef(%TRANS_RXN); 
undef(%RXN_ALTS);  undef(%RXN_DIR);   undef(%RXN_EC); undef(%UPID_RXN); undef(%RXN_UPID);
##################################################################################
#################  	 DONE INPUT PUBCHEM COMPOUNDS		##################
##################################################################################











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
sub CleanNames{
@GREEKS = ("α","β","γ","δ","ε","ζ","η","θ","ι","κ","λ","μ","ν","ξ","ο","π","ρ","ς","σ","τ","υ","φ","χ","ψ","ω");
@GREEKL = ("ALPHA","BETA","GAMMA","DELTA","EPSILON","ZETA","ETA","THETA","IOTA","KAPPA","LAMBDA","MU","NU","XI","OMICRON","PI","RHO","SIGMA","SIGMA","TAU","UPSILON","PHI","CHI","PSI","OMEGA");
        $nameX = $_[0];
        #remove junk punctuation/standardize
        $sta=0; $end=1;
        while($end ne $sta){
                $sta=$nameX;
                #swap greek symbols for text
                for my $g (0..$#GREEKL){ #fix pathbank and other greek symbols
                        if($nameX =~/($GREEKS[$g])/){
                                $nameX =~ s/$GREEKS[$g]/$GREEKL[$g]/g;
                }       }
                $nameX =~ s/\%2B(\d*)/$1\+/g;   #fix html +/- code (HMDB db)
                $nameX =~ s/\%2D(\d*)/$1\-/g;   #fix html +/- code (HMDB db)
                $nameX =~ s/(ARROW|STEREO|RIGHT|LEFT|\-)*\&/\&/g; #fix html +/- code (rhea)
                $nameX =~ s/\&\w+\;\/*//g; #fix html +/- code (rhea)

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
                $nameX =~ s/^(LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|UNCHARACTERIZED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)//g;
                $nameX =~ s/^(UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)//g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|UNCHARACTERIZED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](UNDESCRIBED|UNKNOWN|UNCULTIVATED|UNCULTURED|UNIDENTIFIED|UNCLASSIFIED|UNASSIGNED)[\b\_]/\_/g;
                $nameX =~ s/[\b\_](LIKE|CANDIDATUS|CANDIDATUAS|VOUCHERED|UNCHARACTERIZED|ASSOCIATED|CONTAMINATION|U*N*SCREENED|COMBINED|PUTATIVE)$//g;
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
		$name=CleanNames($name);
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



