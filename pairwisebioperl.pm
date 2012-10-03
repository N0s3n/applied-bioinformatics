package pairwiseAlignment;
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

sub pa {
    my $fhdelta;
    #Temporary filenames for subset of assemblyfile 
    my $subQ_tempy; 
    my $subR_tempy;
    #Temporary filename for deltafile
    my $nucmer_tempy;
    #contigs array, in first level two elements. Element 0 contains all the contigs matching from query sequence.
    #Element 1 contains all matches from reference sequence. The matches are represented as three dimensional arrays containing contig id, start of match and end of match.
    my @contigs;
    my @contigs_pos;
    my @contigs_p;
    my @filenames = @_;
    my @pre_comb;

    for my $qi (0..$#filenames) {
	for my $ri (0..$#filenames) {
	
	    if ($qi != $ri) {
		$nucmer_tempy = "nucmer_out_temp".$qi."_".$ri; 

		qx(nucmer --prefix="$nucmer_tempy" "$filenames[$qi]" "$filenames[$ri]");
		#TODO: error handling if nucmer fails

		#Create the contigs array
		open(my $fhdelta, "<",$nucmer_tempy.".delta") 
		    or die "cannot open < ".$nucmer_tempy.".delta : $!";
		@contigs = get_contigs($fhdelta);
		close($fhdelta);

		#Check positions.
		#if (not(($qi==0)and($ri==1))) {
		#    if ($pre_comb[0]==$qi) {
	#		my @a = position ($contigs_p[0],$contigs[0]);
#			$contigs[0] = [ @a ];		    
#		    }
#		    elsif ($pre_comb[0]==$ri) {
#			my @b = position ($contigs_p[0],$contigs[1]);
#			$contigs[1] = [ @b ];
#			return @b;
#		    }
#
#		    if ($pre_comb[1]==$qi) {
#			my @c = position ($contigs_p[1],$contigs[0]);
#			$contigs[0] = [ @c ];
#		    }
#		    elsif ($pre_comb[1]==$qi) {
#			my @d = position ($contigs_p[1],$contigs[1]);
#			$contigs[1] = [ @d ];
#		    }
#		}
		
		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "sub_tempR".$ri.$qi.".fasta"; 
		#Output now occurs in filter().
		filter($filenames[$qi],$subQ_tempy,@{$contigs[0]});
		filter($filenames[$ri],$subR_tempy,@{$contigs[1]});
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;

		@pre_comb = ($qi,$ri);
		@contigs_p = @contigs;#Store information for position checks in next iteration
	    }
	}
    }

    #qx(rm -f "$nucmer_tempy.delta");
    return @contigs;
}

#function that goes through the .delta file and stores the contig ids and their positions in an array.
sub get_contigs {
    my @contigs;
    my $c1;
    my $c2;

    for my $fh (@_) {
	my $i = 0;
	for my $line (<$fh>) {
	    if ($line =~ /^>([^\s]+)\s+([^\s]+)/) { #Matches the contig ids
		$c1 = $1;
		$c2 = $2;
	    }
	    if ($line =~ /^(\d+)\s(\d+)\s(\d+)\s(\d+)\s\d+\s\d+\s\d+$/) { #Matching of position information
		#Storing contig ids and position information in the contig array
		if ($1 < $2) {
			@{$contigs[0][$i]} = ($c1,$1,$2);
		}
                else {
                        @{$contigs[0][$i]} = ($c2,$2,$1);
                }
		if ($3 < $4) {
			@{$contigs[1][$i]} = ($c2,$3,$4);
		}
		else {
			@{$contigs[1][$i]} = ($c2,$4,$3);
		}
		$i = $i + 1;
	    }
	}
    }
    return @contigs;
}

sub filter {
	my ($filename, $output, @contigs) = @_;
	my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
	my $filtered = Bio::SeqIO->new(-file => ">$output", -format => 'fasta' );
	while (my $seq = $fh->next_seq)
	{
		my $contigname = $seq->display_id;
		my $hits = 0;
		for my $index (0 .. $#contigs)
		{
			if ($contigs[$index][0] eq $contigname)
			{
				$hits++;
				#print $contigs[$index][0]."\t".$contigname."\n";
				#print $index . "\t" . $seq->length . "\t$contigs[$index][2]\n";
				my $subseq = $seq->subseq($contigs[$index][1],$contigs[$index][2]);
				$subseq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname"."\.$hits");
				$filtered->write_seq($subseq);
			}
		}
	}
}
1;
