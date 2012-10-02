package pairwiseAlignment;
use strict;
use warnings;
use lib do {
    use Cwd qw/realpath/;
    realpath($0) =~ /^(.*)\//;
    "$1/BioPerl-1.6.901/install/share/perl";
};
use Bio::Perl;

sub min {
    return (sort {$a<=>$b} @_)[0]; 
}

sub max {
    return (sort {$b<=>$a} @_)[0]; 
}

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
		if (not(($qi==0)and($ri==1))) {
		    if ($pre_comb[0]==$qi) {
			my @a = position ($contigs_p[0],$contigs[0]);
			$contigs[0] = [ @a ];		    
		    }
		    elsif ($pre_comb[0]==$ri) {
			my @b = position ($contigs_p[0],$contigs[1]);
			$contigs[1] = [ @b ];
		    }

		    if ($pre_comb[1]==$qi) {
			my @c = position ($contigs_p[1],$contigs[0]);
			$contigs[0] = [ @c ];
		    }
		    elsif ($pre_comb[1]==$qi) {
			my @d = position ($contigs_p[1],$contigs[1]);
			$contigs[1] = [ @d ];
		    }
		}
		
		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "sub_tempR".$ri.$qi.".fasta"; 
		open(my $fhq, ">",$subQ_tempy) 
		    or die "cannot open < ".$subQ_tempy.": $!";
		open(my $fhr, ">",$subR_tempy) 
		    or die "cannot open < ".$subR_tempy.": $!";
		print $fhq filter($filenames[$qi],@{$contigs[0]});
		print $fhr filter($filenames[$ri],@{$contigs[1]});
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;
		close($fhq);
		close($fhr);

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
		@{$contigs[0][$i]} = ($c1,$1,$2);
		@{$contigs[1][$i]} = ($c2,$3,$4);
		$i = $i + 1;
	    }
	}
    }

    return @contigs;
}

sub filter {
	my ($filename, @contigs) = @_;
	my @keys;
	my $boolean = 0;
	for my $index (0 .. $#contigs)
	{
	    push (@keys, $contigs[$index][0]);
	}
	
	my $found = 0;
	my $content = "";
	open (my $in, "<", "$filename");
	while (my $in = <$in>)
	{
		if($in =~ /^>([^\s]+)/)
		{
			my $contigname = $1;
			$found = 0;
			foreach (@keys)
			{
				if ($_ eq $contigname)
				{
					$found = 1;
					$content .= $in;
					last();
				}
			}
		}
		elsif ($found == 1)
		{
			$content .= $in;
		}
	}
	close($in);
	return $content;
}

sub position {
    my @pos;
    my @c_p = @{ $_[0] };#contig information from previous iteration
    my @c   = @{ $_[1] };#contig information from current iteration

    for my $pre (@c_p) {

	for my $cur (@c) {

	    if ($pre->[0] eq $cur->[0]) {

		my $startP = min($pre->[1],$pre->[2]);
		my $endP = max($pre->[1],$pre->[2]);
		my $rangeP = Bio::Range->new(-start => $startP, -end => $endP);
		my $startC = min($cur->[1],$cur->[2]);
		my $endC = max($cur->[1],$cur->[2]);
		my $rangeC = Bio::Range->new(-start => $startC, -end => $endC);
		my ($start,$end,$strand) = $rangeC->intersection($rangeP);

		if(defined($start)) {
		    push (@pos ,[ $pre->[0],$start,$end ]);
		}

	    }

	}

    }

    return @pos;
}
1;
