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
    my %contigs;
    my @filenames = @_;
    mkdir "pabl";

    for my $qi (0..$#filenames) {
	for my $ri (0..$#filenames) {
	
	    if ($qi != $ri) {
		$nucmer_tempy = "pabl/nucmer_out_temp".$qi."_".$ri; 

		qx(nucmer -l 10 --prefix="$nucmer_tempy.pre" "$filenames[$qi]" "$filenames[$ri]");
		qx(delta-filter -i 99.9 "$nucmer_tempy.pre.delta" > "$nucmer_tempy.delta");
		qx(rm -f $nucmer_tempy.pre.delta);
		#TODO: error handling if nucmer fails

		#Create the contigs array
		open(my $fhdelta, "<",$nucmer_tempy.".delta") 
		    or die "cannot open < ".$nucmer_tempy.".delta : $!";
		%contigs = get_contigs($fhdelta);
		close($fhdelta);
		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "pabl/sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "pabl/sub_tempR".$qi.$ri.".fasta"; 
		#Output now occurs in filter().
		filter($filenames[$qi],$subQ_tempy,%contigs);
		filter($filenames[$ri],$subR_tempy,%contigs);
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;
	    }
	}
    }
    qx(cp $subQ_tempy output/);
    #qx(rm -f "$nucmer_tempy.delta");
    #return %contigs;
}

#function that goes through the .delta file and stores the contig ids and their positions in an array.
sub get_contigs {
    my %contigs;
    my $c1;
    my $c2;

    for my $fh (@_) {
	for my $line (<$fh>) {
	    my @pos;
	    if ($line =~ /^>([^\s]+)\s+([^\s]+)/) { #Matches the contig ids
		$c1 = $1;
		$c2 = $2;
	    }
	    if ($line =~ /^(\d+)\s(\d+)\s(\d+)\s(\d+)\s\d+\s\d+\s\d+$/) { #Matching of position information
		#Storing contig ids and position information in the contig array
		if ($contigs{$c1}) {
                    @pos = @{ $contigs{$c1} };
                }
		if ($1 < $2) {
                	push (@pos,$1,$2);
		}
		else {
			push (@pos,$2,$1);
		}
                $contigs{$c1} = [ @pos ];
		@pos = ();
                if ($contigs{$c2}) {
                    @pos = @{ $contigs{$c2} };
                }
                if ($3 < $4) {
                        push (@pos,$3,$4);
                }
                else {
                        push (@pos,$4,$3);
                } 
                $contigs{$c2} = [ @pos ];
	    }
	}
    }
    return %contigs;
}

sub filter {
	my ($filename, $output, %contigs) = @_;
        #Reads the sequence being filtered:
	my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
        #Opens a new file to write to:
	my $filtered = Bio::SeqIO->new(-file => ">$output", -format => "fasta" );
	#Goes through the contigs one by one:
	while (my $seq = $fh->next_seq)
	{
		my $contigname = $seq->display_id;
		my $hits = 0;
		#tests to see if there are any matches for the current contig:
		if ($contigs{$contigname})
		{
			my $hittotal = $#{$contigs{$contigname}};
			#Goes through all the matches of each contig.
			for (my $index2 = 0; $index2 < $hittotal; $index2 = $index2 + 2)
			{
				#The hit number is concatenated to the contig name to prevent mismatches of contigs from the same origin.
				$hits++;
				#Extracts the matching intervals:
				my $subseq = $seq->subseq($contigs{$contigname}[$index2],$contigs{$contigname}[$index2+1]);
				#An unsuccesful experiment in storing info in the descriptions.
				#my $descr = 0;
				#if ($seq->desc)
				#{
				#	$descr = $seq->desc;
				#	$descr++;
				#}
				if (length($subseq) > 200) #Filters out hits shorter than the threshold (200, in this case).
				{
					$subseq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname"."\.$hits");
					$filtered->write_seq($subseq);
				}
			}
		}
	}
}
1
