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
$kegg_cpds = `wget -q -O - https://rest.kegg.jp/list/compound`;
$kegg_rxns = `wget -q -O - https://rest.kegg.jp/list/reaction`;


#OPEN/CHECK INPUTS:
$time=localtime;
print "OPEN / CHECK INPUTS time $time\n";
open(INTRCH,	"getSubstrates.py")		||die "unable to open getSubstrates.py: $!\n";

#OPEN/CHECK OUTPUTS:
open(OUTKGC, ">", "KEGG_CPD_DB.txt")||die;
open(OUTKGR, ">", "KEGG_RXN_DB.txt")||die;
################################################
################################################

 







##################################################################################
######################   	LOAD UNIVERSAL INFO	 #########################
##################################################################################


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
#######################	  	   INPUT KEGG COMPOUNDS 	##################
##################################################################################
$time=localtime;
print "INPUT KEGG CPDS time $time\n";
@KG_CPD = split("\n", $kegg_cpds);
$ckc=@KG_CPD;
$cnt=0;
foreach my $x (@KG_CPD){
        $cnt++;
        $time=localtime;
        if($cnt%1000==0){print "on $cnt of $ckc time $time\n";}
        (my $cpid, my $names)=split("\t",$x);
        if($cpid !~ /C\d+/){next;}
	$cpid =~ /(C\d+)/; $cpd=$1;
	$CMPD_ALTS{$cpd}{$cpd}="KG";
        $cpd_info = `wget -q -O - https://rest.kegg.jp/get/$cpid`;
	$cpd_info = uc($cpd_info);
	@lines=split("\n",$cpd_info);
	$names=uc($names);
	@names=split(";",$names);
	foreach my $name (@names){ $name=CleanNames($name); 	 $CMPD_NAME{$cpd}{$name}="KG"; } 
	foreach my $i (@lines){
		if($i=~/^FORMULA\s+(\S+)\;*/){	 $form=$1; 	 $CMPD_FORM{$cpd}=CleanNames($form);}
		if($i=~/^EXACT_MASS\s+(\S+)\;*/){$mass=$1; 	 $CMPD_MASS{$cpd}=$mass;}
		if($i=~/^\s+PUBCHEM\:\s+(\d+)/){ $alt="CID:".$1; $CMPD_ALTS{$cpd}{$alt}="KG"; }
		if($i=~/^\s+CHEBI\:\s+(\d+)/){ $alt="CHEBI:".$1; $CMPD_ALTS{$cpd}{$alt}="KG"; }
	}
	#leaving the hmdb id code below - takes 4X as long so not using it ATM
	#$hmid = `wget -q -O - https://www.genome.jp/dbget-bin/get_linkdb?-t+hmdb+$cpid`;
	#if($hmid =~ /(HMDB\d+)/){ $alt=$1; while(length($alt)<11){$alt=~s/HMDB/HMDB0/;} $CMPD_ALTS{$cpd}{$alt}="KG"; }
}




#OUTPUT KEGG CPDS
$time=localtime;
print "OUTPUT KEGG CPDS time $time\n";
print OUTKGC "cpd\tformula\tmass\tcharge\tdb_src\ttcdbs\tnames\tkeggcpd\tchebcpd\thmdbcpd\tpubccpd\tinchcpd\tbioccpd\n";
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
	if(keys %TCDBS > 0){ foreach my $tc (sort(keys %TCDBS)){ push(@TCDBS,$tc); } @TCDBS=nsort(@TCDBS); $tcdb=join(";",@TCDBS);}}
	$hmdbcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $hmdbcpd); @CPDS=nsort(@CPDS); $hmdbcpd=join(";",@CPDS);
        $pubccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $pubccpd); @CPDS=nsort(@CPDS); $pubccpd=join(";",@CPDS);
	$inchcpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $inchcpd); @CPDS=nsort(@CPDS); $inchcpd=join(";",@CPDS);
	$bioccpd=~s/\;+$//g; @CPDS=(); @CPDS=split(";", $bioccpd); @CPDS=nsort(@CPDS); $bioccpd=join(";",@CPDS);
	print OUTKGC "$cpd\t$form\t$mass\t$char\t$sorc\t$tcdb\t$name\t$keggcpd\t$chebcpd\t$hmdbcpd\t$pubccpd\t$inchcpd\t$bioccpd\n";
}
#####################################################################################
#######################	         DONE INPUT KEGG COMPOUNDS	#####################
#####################################################################################



##################################################################################
#######################	  	   INPUT KEGG REACTIONS 	##################
##################################################################################
$time=localtime;
print "INPUT KEGG RXNS time $time\n";
@KG_RXN = split("\n", $kegg_rxns);
$rkc=@KG_RXN;
$cnt=0;
foreach my $x (@KG_RXN){
	$cnt++;
	$time=localtime;
	if($cnt%1000==0){print "on $cnt of $rkc time $time\n";}
        @stuff=split("\t",$x);
	$rxid=$stuff[0];
        if($rxid !~ /R\d+/){next;}
	$rxn_info = `wget -q -O - https://rest.kegg.jp/get/$rxid`;
	$rxn_info = uc($rxn_info);
	@LEFTS = ();
	@RITES = ();
	$rxn=''; $ec=''; $l=''; $r=''; $rhea='';
	@lines=split("\n",$rxn_info);
	foreach my $i (@lines){
		if($i=~/^ENTRY\s+(R\d+)/){$rxn=$1;}
		if($i=~/^EQUATION\s+(.*)\=(.*)$/){
			$l=$1; $r=$2;
			@LEFTS = ( $l=~/(C\d+)/g );
			@RITES = ( $r=~/(C\d+)/g );
		}
		if($i=~/^ENZYME\s+(\d+\.\d+\.\d+\.N*\d+)/){
			@ES = ($i=~/(\d+\.\d+\.\d+\.N*\d+)/g);    #keep the first one with UPIDS
			foreach my $ex (@ES){ if(exists($EC2UPID{$ex})){$RXN_EC{$rxn}{$ex}="KG"; $ec=$ex; last;}}
			if($ec !~/\w/){$ec=$ECS[0];}
		}
		if($i=~/^DBLINKS\s+RHEA\:\s*(\d+)/){$rhea=$1; }
	}
	if($rxn !~/R\d+/){next;}
	$RXN_ALTS{$rxn}{$rxn}="KG";
	$RXN_ALTS{$rxn}{$rhea}="KG";
        $RXN_DIR{$rxn}="LTR";

	#get unprot from ec:
	if(keys %{$EC2UPID{$ec}} < 1 && $ec=~/\d/){ #enzyme db didn't have upids, maybe kegg
		@UPID=();
		$ecup_info=`wget -q -O - https://www.genome.jp/dbget-bin/get_linkdb?-t+uniprot+ec:$ec`;
		if($ecup_info=~/No link information was found/i){next;}
		@UPID=($ecup_info =~ /entry.up.(\w+)\"/ig);
		foreach my $upid (@UPID){ 
			$RXN_UPID{$rxn}{$upid}="KG"; 
			$UPID_RXN{$upid}{$rxn}="KG"; 
		}
	}
	elsif(keys %{$EC2UPID{$ec}} >= 1 && $ec=~/\d/){
		foreach my $upid (keys %{$EC2UPID{$ec}}){ 
		$RXN_UPID{$rxn}{$upid}="KG";
		$UPID_RXN{$upid}{$rxn}="KG";}}
	else{}
	foreach my $left (@LEFTS){ if($left!~/\w/){next;}
		$CMPD_ALTS{$left}{$left}="KG"; 
		$LRXN_CPD{$rxn}{$left}="INSIDE";  
		$LCPD_RXN{$left}{$rxn}="INSIDE";}
	foreach my $right (@RITES){if($right!~/\w/){next;} 
		$CMPD_ALTS{$right}{$right}="KG";
		$RRXN_CPD{$rxn}{$right}="INSIDE"; 
		$RCPD_RXN{$right}{$rxn}="INSIDE";}
}


#######################################################	
#OUTPUT KEGG RXN DATA
#######################################################	
$cnt=0;
$time=localtime;
print "OUTPUT Kegg rxn time $time\n";
print OUTKGR "rxn\tdb_src\trxn_dir\talt_rhea\talt_kegg\talt_bioc\t";
print OUTKGR "rhea_lcpd\trhea_lloc\trhea_ltrn\trhea_rcpd\trhea_rloc\trhea_rtrn\t";
print OUTKGR "kegg_lcpd\tkegg_lloc\tkegg_ltrn\tkegg_rcpd\tkegg_rloc\tkegg_rtrn\t";
print OUTKGR "bioc_lcpd\tbioc_lloc\tbioc_ltrn\tbioc_rcpd\tbioc_rloc\tbioc_rtrn\tUP_ID\tEC\n";
foreach my $rxn (sort(keys %RXN_ALTS)){
       	@rhea_lcpd=();  @kegg_lcpd=();  @bioc_lcpd=();
        @rhea_lloc=();  @kegg_lloc=();  @bioc_lloc=();
        @rhea_ltrn=();  @kegg_ltrn=();  @bioc_ltrn=();
        @rhea_rcpd=();  @kegg_rcpd=();  @bioc_rcpd=();
        @rhea_rloc=();  @kegg_rloc=();  @bioc_rloc=();
        @rhea_rtrn=();  @kegg_rtrn=();  @bioc_rtrn=();

        $rhea_lcpd='';  $kegg_lcpd='';  $bioc_lcpd='';
        $rhea_lloc='';  $kegg_lloc='';  $bioc_lloc='';
        $rhea_ltrn='';  $kegg_ltrn='';  $bioc_ltrn='';
        $rhea_rcpd='';  $kegg_rcpd='';  $bioc_rcpd='';
        $rhea_rloc='';  $kegg_rloc='';  $bioc_rloc='';
        $rhea_rtrn='';  $kegg_rtrn='';  $bioc_rtrn='';
        $rrx='';        $krx='';        $brx='';
        $dir=$RXN_DIR{$rxn};
        $sorc=$RXN_ALTS{$rxn}{$rxn};


        foreach my $cpd (sort(keys %{$LRXN_CPD{$rxn}})){
                if($cpd !~/\w/){next;}
                if($TRANS_RXN{$rxn}{$cpd}=~/PORT/){ $tdir=$TRANS_RXN{$rxn}{$cpd};} else{$tdir="NOPORT";} #import/export/biport/noport
                if( $LRXN_CPD{$rxn}{$cpd}=~/SIDE/){ $tloc=$LRXN_CPD{$rxn}{$cpd};}  else{$tloc="INSIDE";} #inside or outside
                foreach my $alt (sort(keys %{$CMPD_ALTS{$cpd}})){
                  #put in a clean namess sub for cpds - also why duplicate compounds? with hash how?
                        print "rxn $rxn cpd $cpd alt $alt\n";
                           if($alt=~/^INCHI\:\w+|^HMDB\d+|^CID\:\d+/){next;}
                        elsif($alt=~/^C\d+$/){          push(@kegg_lcpd,$alt); push(@kegg_lloc,$tloc); push(@kegg_ltrn,$tdir);}
                        elsif($alt=~/^CHEBI.\d+$/){     push(@rhea_lcpd,$alt); push(@rhea_lloc,$tloc); push(@rhea_ltrn,$tdir);}
                        else{ $alt=CleanNames($alt);    push(@bioc_lcpd,$alt); push(@bioc_lloc,$tloc); push(@bioc_ltrn,$tdir);}
                }
        }
        foreach my $cpd (sort(keys %{$RRXN_CPD{$rxn}})){
                if($cpd !~/\w/){next;}
                if($TRANS_RXN{$rxn}{$cpd}=~/PORT/){ $tdir=$TRANS_RXN{$rxn}{$cpd};} else{$tdir="NOPORT";} #import/export/biport/noport
                if( $RRXN_CPD{$rxn}{$cpd}=~/SIDE/){ $tloc=$RRXN_CPD{$rxn}{$cpd}; } else{$tloc="INSIDE";} #inside or outside
                foreach my $alt (sort(keys %{$CMPD_ALTS{$cpd}})){
                           if($alt=~/^INCHI\:\w+|^HMDB\d+|^CID\:\d+/){next;}
                        elsif($alt=~/^C\d+$/){          push(@kegg_rcpd,$alt); push(@kegg_rloc,$tloc); push(@kegg_rtrn,$tdir);}
                        elsif($alt=~/^CHEBI.\d+$/){     push(@rhea_rcpd,$alt); push(@rhea_rloc,$tloc); push(@rhea_rtrn,$tdir);}
                        else{ $alt=CleanNames($alt);    push(@bioc_rcpd,$alt); push(@bioc_rloc,$tloc); push(@bioc_rtrn,$tdir);}

                }
        }

        $rhea_lcpd=join(";",@rhea_lcpd);  $kegg_lcpd=join(";",@kegg_lcpd);  $bioc_lcpd=join(";",@bioc_lcpd);
        $rhea_lloc=join(";",@rhea_lloc);  $kegg_lloc=join(";",@kegg_lloc);  $bioc_lloc=join(";",@bioc_lloc);
        $rhea_ltrn=join(";",@rhea_ltrn);  $kegg_ltrn=join(";",@kegg_ltrn);  $bioc_ltrn=join(";",@bioc_ltrn);
        $rhea_rcpd=join(";",@rhea_rcpd);  $kegg_rcpd=join(";",@kegg_rcpd);  $bioc_rcpd=join(";",@bioc_rcpd);
        $rhea_rloc=join(";",@rhea_rloc);  $kegg_rloc=join(";",@kegg_rloc);  $bioc_rloc=join(";",@bioc_rloc);
        $rhea_rtrn=join(";",@rhea_rtrn);  $kegg_rtrn=join(";",@kegg_rtrn);  $bioc_rtrn=join(";",@bioc_rtrn);

        %ECS=(); %UPIDS=(); @UPIDS=(); $ec=''; $ex='';
        foreach my $ec   (keys %{$RXN_EC{$rxn}}){       $ECS{$ec}++;
          foreach my $upid (keys %{$EC2UPID{$ec}}){     $UPIDS{$upid}++;}}
        foreach my $upid (keys %{$RXN_UPID{$rxn}}){     $UPIDS{$upid}++;}
        foreach my $ex (sort(keys %ECS)){ if($ex=~/\d/){$ec = $ex; last;} }
        foreach my $upid (keys %UPIDS){ push(@UPIDS,$upid); }
        @UPIDS=nsort(@UPIDS);
        $upid=join(";",@UPIDS);

        foreach my $alt (sort(keys %{$RXN_ALTS{$rxn}})){
                   if($alt =~ /^\d+$/){ $rrx.=$alt.";"; }
                elsif($alt =~ /^R\d+$/){$krx.=$alt.";"; }
                elsif($alt =~ /RXN/){   $brx.=$alt.";"; }
                else{}
        }
        $one="$rxn\t$sorc\t$dir\t$rrx\t$krx\t$brx\t";
        $two="$rhea_lcpd\t$rhea_lloc\t$rhea_ltrn\t$rhea_rcpd\t$rhea_rloc\t$rhea_rtrn\t";
        $thr="$kegg_lcpd\t$kegg_lloc\t$kegg_ltrn\t$kegg_rcpd\t$kegg_rloc\t$kegg_rtrn\t";
        $for="$bioc_lcpd\t$bioc_lloc\t$bioc_ltrn\t$bioc_rcpd\t$bioc_rloc\t$bioc_rtrn\t$upid\t$ec";
        $out="$one$two$thr$for";
        $out=~s/\;+\t/\t/g;
        $out=~s/\;+$//g;


	print OUTKGR "$out\n";
}
undef(%CMPD_ALTS); undef(%CMPD_CHAR); undef(%CMPD_FORM); undef(%CMPD_MASS); undef(%CMPD_NAME); 
undef(%LRXN_CPD);  undef(%RRXN_CPD);  undef(%TRANS_CPD); undef(%LCPD_RXN);  undef(%RCPD_RXN);  undef(%TRANS_RXN); 
undef(%RXN_ALTS);  undef(%RXN_DIR);   undef(%RXN_EC); undef(%UPID_RXN); undef(%RXN_UPID);
###################################################################################
##############     DONE INPUT KEGG  REACTIONS   ####################
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
return($nameX);}


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



