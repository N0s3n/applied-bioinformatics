package pairwisebioperl;
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
		$nucmer_tempy = "output/nucmer_out_temp".$qi."_".$ri; 

		qx(nucmer -l 10 --prefix="$nucmer_tempy.pre" "$filenames[$qi]" "$filenames[$ri]");
		qx(delta-filter -i 99.9 "$nucmer_tempy.pre.delta" > "$nucmer_tempy.delta");
		qx(rm -f $nucmer_tempy.pre.delta);
		#TODO: error handling if nucmer fails

		#Create the contigs array
		open(my $fhdelta, "<",$nucmer_tempy.".delta") 
		    or die "cannot open < ".$nucmer_tempy.".delta : $!";
		@contigs = get_contigs($fhdelta);
		close($fhdelta);

		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "output/sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "output/sub_tempR".$qi.$ri.".fasta"; 
		#Output now occurs in filter().
		filter($filenames[$qi],$subQ_tempy,@{$contigs[0]});
		filter($filenames[$ri],$subR_tempy,@{$contigs[1]});
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;
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
                        @{$contigs[0][$i]} = ($c1,$2,$1);
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
	#Reads the sequence being filtered:
	my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
	#Opens a new file to write to:
	my $filtered = Bio::SeqIO->new(-file => ">$output", -format => "fasta" );
	while (my $seq = $fh->next_seq)
	{
		my $contigname = $seq->display_id;
		my $hits = 0;
		for my $index (0 .. $#contigs)
		{
			#Tests the contig names in the delta files versus the contig names in the sequence being filtered:
			if ($contigs[$index][0] eq $contigname)
			{
				#The three lines below exist for testing purposes:
				#print $#contigs . "\t";
				#print $contigs[$index][0]."\t".$contigname."\n";
				#print $index . "\t" . $seq->length . "\t$contigs[$index][2]\n"; 
				#The hit number is concatenated to the contig name.
				$hits++;
				#Extracts the subsequences that match:
				my $subseq = $seq->subseq($contigs[$index][1],$contigs[$index][2]);
				#An unsuccesful experiment in storing info in the descriptions.
				#my $descr = 0;
				#if ($seq->desc)
				#{
				#	$descr = $seq->desc;
				#	$descr++;
				#}
				if (length($subseq) > 200) #Filters out hits shorter than the threshold (200, in this case).
				{
					#Writes a new sequence to the $output file:
					$subseq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname"."\.$hits");
					$filtered->write_seq($subseq);
				}
			}
		}
	}
}
1;
