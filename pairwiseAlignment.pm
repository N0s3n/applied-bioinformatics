package pairwiseAlignment;
use strict;
use warnings;
use lib do {
    use Cwd qw/realpath/;
    realpath($0) =~ /^(.*)\//;
    "$1/BioPerl-1.6.901/install/share/perl";
};
use Bio::Perl;
use Data::Dumper;
use Bio::SeqIO;

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
    my (%contigs1,%contigs2);
    my @contigs_pos;
    my (%contigs_p1,%contigs_p2);
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
		my ($c1,$c2) = get_contigs($fhdelta);
		%contigs1 = %$c1;
		%contigs2 = %$c2;
		close($fhdelta);

		#Check positions.
		if (not(($qi==0)and($ri==1))) {
		    if ($pre_comb[0]==$qi) {
			#my %a = position (\%contigs_p1,\%contigs1);
			%contigs1 =  position (\%contigs_p1,\%contigs1);		    
		    }
		    elsif ($pre_comb[0]==$ri) {
			#my %b = position (\%contigs_p1,\%contigs2);
			%contigs2 = position (\%contigs_p1,\%contigs2);;
		    }

		    if ($pre_comb[1]==$qi) {
			#my %c = position (\%contigs_p2,\%contigs1);
			%contigs1 = position (\%contigs_p2,\%contigs1);;
		    }
		    elsif ($pre_comb[1]==$qi) {
			#y %d = position (\%contigs_p2,\%contigs2);
			%contigs2 = position (\%contigs_p2,\%contigs2);;
		    }
		}
		
		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "sub_tempR".$ri.$qi.".fasta"; 
		# open(my $fhq, ">",$subQ_tempy) 
		#     or die "cannot open < ".$subQ_tempy.": $!";
		# open(my $fhr, ">",$subR_tempy) 
		#     or die "cannot open < ".$subR_tempy.": $!";
#		print $fhq 
		filter($filenames[$qi],$subQ_tempy,\%contigs1);
#		print $fhr 
		filter($filenames[$ri],$subR_tempy,\%contigs2);
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;
		# close($fhq);
		# close($fhr);

   open(my $fhtest, ">>","test_contig_array".$qi."_".$ri.".txt") 
	or die "cannot open < ".$nucmer_tempy.".delta : $!";
		print $fhtest "Contig:\n";
		print $fhtest Dumper \%contigs1;
		print $fhtest "\n";
		close($fhtest);

		@pre_comb = ($qi,$ri);
		%contigs_p1 = %contigs1;#Store information for position checks in next iteration
		%contigs_p2 = %contigs2;#Store information for position checks in next iteration
	    }
	}
    }

    #qx(rm -f "$nucmer_tempy.delta");
    return \%contigs1,\%contigs2;
}

#function that goes through the .delta file and stores the contig ids and their positions in an array.
sub get_contigs {
    #my @contigs;
    my (%contigs1,%contigs2);
    my ($c1,$c2);

    for my $fh (@_) {

	for my $line (<$fh>) {

	    my (@pos1,@pos2);
	    if ($line =~ /^>([^\s]+)\s+([^\s]+)/) { #Matches the contig ids
		$c1 = $1;
		$c2 = $2;
	    }
	    if ($line =~ /^(\d+)\s(\d+)\s(\d+)\s(\d+)\s\d+\s\d+\s\d+$/) { #Matching of position information
		#Storing contig ids and position information in the contig array

		if (defined($contigs2{$c1})) {
		    @pos1 = @{ $contigs1{$c1} };
		}
		if (defined($contigs2{$c2})) {
		    @pos2 = @{ $contigs2{$c2} };
		}

		push (@pos1, [$1,$2]);
		$contigs1{$c1} = [ @pos1 ];
		push (@pos2,[$3,$4]);
		$contigs2{$c2} = [ @pos2 ];
	    }
	}
    }

    return \%contigs1,\%contigs2;
}

sub filter {
	my ($filename, $output, $cont) = @_;
	my @contigs = keys (%{$cont});
	my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
	my $filtered = Bio::SeqIO->new(-file => ">$output", -format => 'fasta' );
	while (my $seq = $fh->next_seq)
	{
		my $contigname = $seq->display_id;
		my $hits = 0;
		for my $cid (@contigs)
		{
			if ($cid eq $contigname)
			{
				$hits++;
				#my $subseq = $seq->subseq($contigs[$index][1],$contigs[$index][2]);
				#$subseq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname"."\.$hits");
				$filtered->write_seq($seq);
			}
		}
	}
}

sub position {
    my %cpos;
    my %c_p = %{ $_[0] };#contig information from previous iteration
    my %c   = %{ $_[1] };#contig information from current iteration

    for my $cid (keys(%c)) {
	for my $cidp (keys(%c_p)) {
	    if (defined($c{$cid})) {
		my @pos;
		for my $po_cur (@{ $c{$cid} }) {
		    if (defined($c_p{$cid})) {
			for my $po_pre (@{ $c_p{$cidp} }) {
			    my $startP = min($po_pre->[0],$po_pre->[1]);
			    my $endP = max($po_pre->[0],$po_pre->[1]);
			    my $rangeP = Bio::Range->new(-start => $startP, -end => $endP);
			    my $startC = min($po_cur->[0],$po_cur->[1]);
			    my $endC = max($po_cur->[0],$po_cur->[1]);
			    my $rangeC = Bio::Range->new(-start => $startC, -end => $endC);
			    my ($start,$end,$strand) = $rangeC->intersection($rangeP);

			    if(defined($start)) {
				push (@pos ,[ $start,$end ]);
				$cpos{$cid} = [ @pos ];
			    }
			}
			
		    }
		}
	    }
	}
    }

    # for my $pre (@c_p) {

    # 	for my $cur (@c) {

    # 	    if ($pre->[0] eq $cur->[0]) {

    # 		my $startP = min($pre->[1],$pre->[2]);
    # 		my $endP = max($pre->[1],$pre->[2]);
    # 		my $rangeP = Bio::Range->new(-start => $startP, -end => $endP);
    # 		my $startC = min($cur->[1],$cur->[2]);
    # 		my $endC = max($cur->[1],$cur->[2]);
    # 		my $rangeC = Bio::Range->new(-start => $startC, -end => $endC);
    # 		my ($start,$end,$strand) = $rangeC->intersection($rangeP);

    # 		if(defined($start)) {
    # 		    push (@pos ,[ $pre->[0],$start,$end ]);
    # 		}

    # 	    }

    # 	}

    # }

    return %cpos = %cpos;
}
1;
