use warnings;

$bioc_dir='/geomicro/data2/kiledal/UMRAD/Universal_Biological_Compounds_Database/BIOCYC_NF/*';
#$bioc_dir='/geomicro/data2/tealfurn/URDB/UNIPROT/FUNCTIONS/BIOCYC/*';

$x='';
$x = qx{ls -d $bioc_dir};
@BDIRS =split("\n", $x);


foreach my $bdir (@BDIRS){
        $edatf=$bdir."/enzrxns.dat";
	if(-f $edatf && -s $edatf){push(@ENZ_DAT,$edatf);}
}

###########################################################
###########    INPUT BIOCYC PROTEIN INFO    ##############
print "INPUT BIOCYC enzrxns.dat\n";
open(OUTPUT,">","BIOCYC_MONO_RXNS.txt")||die;
$on=0;
$ac=@ENZ_DAT;
foreach my $inenz (@ENZ_DAT){
	print "on $on of $ac $inenz\n";
	open(INENX, $inenz)||die;
	$inrxn=0;
	$on++;
	while(<INENX>){
		$_ = uc($_);
		$_ =~ s/[\r\n]+//;
		if($_ =~ /UNIQUE-ID/){$inrxn=1; %MONO=(); %RXNS=();}
		if($inrxn==1){
		        if($_ =~ /ENZYME\s\-\s(\S*MONOMER\S*)/){$MONO{$1}=1;}
		        if($_ =~ /REACTION\s\-\s(\S+)/){$RXNS{$1}=1;}
		}
		if($_ =~ /^\/\/\s*$/){
		        $inrxn=0;
		        foreach my $mono (keys %MONO){
		                foreach my $rxn (keys %RXNS){
		                        if($rxn=~/RXN/ && $mono =~/MONOMER/){ print OUTPUT "$mono\t$rxn\n"; }
}	}       }       }       }

##################################################################################################################
##################################################################################################################
##################################################################################################################



