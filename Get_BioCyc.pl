if($ARGV[0]!~/\w/){
        print "Please specify your biocyc licence email for user id, distribution password (not biocyc login password!) and biocyc file list.\n";
        print "ex. perl Get_BioCyc.pl pwd=data-12345 usr=biocyc-flatfiles in=biocyc_files.txt \n";
        print "See GitHub for details: https://github.com/TealFurnholm/Universal_Biological_Compounds_Database/wiki \n";
}

print "GETTING ALL INPUT FILES\n";
qx{mkdir BIOCYC_NF};

$argv = join("\t",@ARGV);
$argv =~ /pwd[\s\=]+(\S+)/;     $pwd=$1;
$argv =~ /usr[\s\=]+(\S+)/;     $usr=$1;
$argv =~ /in[\s\=]+(\S+)/;      $inf=$1;
open(INFILE, $inf)||die "unable to open biocyc files $inf: $!\n";


while(<INFILE>){
        if($_ !~/\w/){next;}
        $_ =~ s/[\r\n]+//;
        $path=$_;
        $file=$_;
        $file=~s/.*\///;
        if($file =~ /(biopax|tier|flatfiles)/i){next;}
        $file =~ /([^\/]+)\.tar\.gz/;
        $fold=$1;
        if($fold !~ /\w/){next;}
        $down=0;
        if(!-s "BIOCYC_NF/$fold/rxn-list"){$down=1;}
        if(!-s "BIOCYC_NF/$fold/enzrxns.dat"){$down=1;}
        if(!-s "BIOCYC_NF/$fold/reactions.dat"){$down=1;}
        if(!-s "BIOCYC_NF/$fold/compounds.dat"){$down=1;}
        #if(!-s "BIOCYC_NF/$fold/protein-links.dat"){$down=1;}
        if(!-s "BIOCYC_NF/$fold/protein-seq-ids.dat"){$down=1;}
        if(!-s "BIOCYC_NF/$fold/metabolic-reactions.xml"){$down=1;}
        if($down==1){
                print "retrieving $file folder $fold path $path\n";
                qx{mkdir BIOCYC_NF/$fold};
                qx{wget -N -q --user=$usr --password=$pwd $path};
                qx{mv $file BIOCYC_NF/$file};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*reactions.dat"};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*compounds.dat"};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*metabolic-reactions.xml"};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*rxn-list"};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*protein-seq-ids.dat"};
                #qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*protein-links.dat"};
                qx{tar -zxf BIOCYC_NF/$file --keep-newer-files --transform='s|.*\/||' --directory="BIOCYC_NF/$fold" --wildcards "*enzrxns.dat"};
                qx{rm "BIOCYC_NF/$file"};
        }
        if($on%10==0){$time=localtime; print "on $on time $time file $file\n";}
        $on++;
}





