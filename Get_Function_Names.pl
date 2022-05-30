use warnings;
$time = localtime;

qx{wget -N https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kog};
qx{wget -N http://current.geneontology.org/ontology/go.obo};
qx{wget -N https://ftp.expasy.org/databases/enzyme/enzyme.dat};
qx{wget -N https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list};
qx{wget -N http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz; gunzip -f Pfam-A.clans.tsv.gz};
qx{wget -N https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/cog-20.def.tab};
qx{wget -N https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv};
$keggrxns  = `wget -q -O - http://rest.kegg.jp/list/reaction`;



open(INGO, "go.obo")||die;
while(<INGO>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	if($_=~/^ID\:\s+(GO.\d+)/){$inrn=1; @IDS=(); push(@IDS,$1); $name='';}
	if($inrn==1 && $_=~/^NAME.\s+(\w.+)/){$name=CleanNames($1);}
	if($inrn==1 && $_=~/^ALT_ID.\s+(GO.\d+)/){push(@IDS,$1);} 
	if($_=~/\[TERM\]/){
		$inrn=0;
		$rc = @IDS;
		if($rc > 0){
			foreach my $rn (@IDS){
				if($name =~ /\w/){$FUNCTION_NAMES{$rn}=$name;}
			}
		}
	}
}

open(INKOG, "kog")||die;
while(<INKOG>){
	if($_!~/^\[\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	if($_=~/^\[([A-Z]+)\]\s+(KOG\d+)\W+(\w.*)/){
		$rn=$2;
		$name=CleanNames($3);
		$name=$1.";".$name;
		if($name =~ /\w/){$FUNCTION_NAMES{$rn}=$name;}
}	}

open(INEC, "enzyme.dat")||die;
while(<INEC>){
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	if($_ =~ /^ID\s+([\d\.\-]+)/){ $rn=$1; $inrn=1; $name='';}
	if($inrn==1 && $_ =~ /^DE\s+(\w.+)/){ $name=CleanNames($1); }
	if($_=~/\/\//){ $inrn=0; 
		if($rn =~ /[\d\.\-]+/ && $name =~ /\w/){$FUNCTION_NAMES{$rn}=$name;}
	}
}

open(INIPR, "entry.list")||die;
while(<INIPR>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	$rn=$stuff[0];
	$name=CleanNames($stuff[2]);
	$FUNCTION_NAMES{$rn}=$name;
}

open(INPFAM, "Pfam-A.clans.tsv")||die;
while(<INPFAM>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	$rn=$stuff[0];
	@stuff=@stuff[2..4];
	%NM=(); @NM=();
	foreach my $n (@stuff){ $nm=CleanNames($n); $NM{$nm}=1; }
	foreach my $n (keys %NM){if($n=~/\w{4}/){push(@NM,$n);}}
	$name=join(";",@NM);
	$FUNCTION_NAMES{$rn}=$name;
}

open(INCOG, "cog-20.def.tab")||die;
while(<INCOG>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	$rn=$stuff[0];
	$name=CleanNames($stuff[2]);
	$name=$stuff[1].";".$name;
	$FUNCTION_NAMES{$rn}=$name;
}

open(INTIGR, "hmm_PGAP.tsv")||die;
while(<INTIGR>){
	if($_!~/\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_);
	if($stuff[1]!~/^TIGR/){next;}
	$rn=$stuff[1];
	$name=CleanNames($stuff[10]);
	$FUNCTION_NAMES{$rn}=$name;
}

@KEGGCD = split("\n", $keggrxns);
%KEGGS=();
foreach my $x (@KEGGCD){
	if($x !~ /\w/){next;}
        $x=uc($x);
	$x=~s/[\r\n]//;
        (my $rn, my $name)=split("\t",$x);
        $rn=~s/RN\://;
	$name=~s/[^\;]+\<\=\>[^\;]+//g;
	@names=split(";",$name);
	for my $i (0..$#names){$names[$i]=CleanNames($names[$i]);}
	$name=join(";",@names);
	if($name !~ /\w/){next;}
	$FUNCTION_NAMES{$rn}=$name;
}

open(OUTPUT, ">", "Function_Names.txt")||die;
foreach my $func (keys %FUNCTION_NAMES){
	$name = $FUNCTION_NAMES{$func};
	print OUTPUT "$func\t$name\n";
}


sub CleanNames{
	$name = $_[0];

	#remove pointless ambiguators
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE|VOUCHERED|UNDESCRIBED|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/UNCHARACTERIZED\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|HYPOTHETICAL)\s/UNCHARACTERIZED\_/g;
	$name =~ s/(UNCHARACTERIZED\_)+/UNCHARACTERIZED\_/g;
	$name =~ s/\-*LIKE\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;

        return($name);
}

