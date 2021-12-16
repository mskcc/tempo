if [ ! -f ccdsGene.txt ] ; then 
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ccdsGene.txt.gz
	gunzip ccdsGene.txt.gz
fi 

process_ccdsgene='
	while (<>) {
		chomp;
        	@F = split /\t/;
        	$F[2] =~ s/^chr//;
        	@exon_starts = split /,/, $F[9];
	        @exon_ends   = split /,/, $F[10];
		$cds_start = $F[6];
		$cds_end = $F[7];
		$strand = $F[3];
                $i = 0;
	        for (0..$#exon_starts) {
			if ($strand == "+") { $type="donor";} else { $type="acceptor"; }
	        	if ($exon_starts[$_] <= $cds_end && $exon_ends[$_] >= $cds_start) {
			print join(
				"\t",
				$F[2],
				$exon_starts[$_] - 3, 
				$exon_starts[$_] + 8,
				"$F[1].$i:$type",
				".",
				"$strand",
			);
			print "\n";
			if ($strand == "+") { $type="acceptor";} else { $type="donor"; }
	                print join(
                        	"\t",
		                $F[2],
	                        $exon_ends[$_] - 8,
			        $exon_ends[$_] + 3,
		                "$F[1].$i:$type",
				".",
				"$strand",
	                );
		        print "\n";
	                $i++;
		}
	}
	}
	'
	perl -e "$process_ccdsgene" ccdsGene.txt > splice_sites.bed

