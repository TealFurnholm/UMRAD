use warnings;


#process
#kegg orgs->genes->kos->rxns
################################################
# GET KEGG ORG GENES AND LINK THRU KOS TO RXNS #
################################################

print "INPUT KEGG KOs ~ RXNs \n";
$ko2rn    = `wget -q -O - http://rest.kegg.jp/link/rn/ko`;
@KO2RN = split("\n", $ko2rn);
foreach my $x (@KO2RN){
        $x=uc($x);
        (my $ko, my $rn)=split("\t",$x);
        $ko =~ s/KO\://;
        $rn =~ s/RN\://;
        if($ko=~/\d/ && $rn=~/\d/){$KO_KRXN{$ko}{$rn}=1;}
}

print "GET KEGG GENE KOS\n";
open(OUTPUT, ">", "KEGG_GENES_RXN.txt")||die "unable to open KEGG_GENES_RXN.txt: $!\n";
$keggorg  = `wget -q -O - http://rest.kegg.jp/list/organism`;
@Y=split("\n",$keggorg); #all kegg organisms
$count=0; 
$totkog=@Y;
foreach my $x (@Y){
        $count++;
        @stuff=split('\s+',$x);
	$org=$stuff[1];
	print "on $count of $totkog org $org\n";
        $z = `wget -q -O - http://rest.kegg.jp/link/ko/$org`;
        $z = uc($z);
        @Z = split("\n",$z);
        foreach my $i (sort(@Z)){
                if($i !~ /K\d\d\d\d\d/){next;}
                $i=~s/\s+KO\:/\t/;
                (my $gene, my $ko)=split("\t",$i);
                $GENE_KO{$gene}{$ko}=1;
		foreach my $rn (sort(keys %{$KO_KRXN{$ko}})){
			if($gene=~/\w/ && $rn =~/\d/){ print OUTPUT "$gene\t$rn\n";}
		}
}       }

$gk=keys %GENE_KO;
$kx=keys %KO_KRXN;
print "gene-ko $gk ko-rxn $kx \n";
###########################################################
###########################################################



