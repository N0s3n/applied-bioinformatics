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

sub savey {
    my ($filename,$thingy) = @_;

    open(my $fhtest, ">>",$filename) 
	or die "cannot open < ".$filename." : $!";
    print $fhtest "Thingy:\n";
    print $fhtest Dumper $thingy;
    print $fhtest "\n";
    close($fhtest);
}

sub pa {
    my $fhdelta;
    #Temporary filenames for subset of assemblyfiles 
    my $subQ_tempy; 
    my $subR_tempy;
    #Temporary filename for deltafile
    my $nucmer_tempy;
    my (%contigs1,%contigs2);
    my @contigs_pos;
    my (%contigs_p1,%contigs_p2);
    my @filenames = @_;
    my @pre_comb;#keeps track of which assemblies were combined in the previous iteration
    my $i = 0;
    my @names;

    for my $qi (0..$#filenames) {
	for my $ri (0..$#filenames) {
	
	    if ($qi != $ri) {
		$nucmer_tempy = "unf_nucmer_out_temp".$qi."_".$ri; 
		$names[$i] = $nucmer_tempy.".delta";

		qx(nucmer --prefix="$nucmer_tempy" "$filenames[$qi]" "$filenames[$ri]");
		$nucmer_tempy = "nucmer_out_temp".$qi."_".$ri.".delta"; 
		qx(delta-filter -i 99 -l 100 $names[$i] > "$nucmer_tempy");		
		qx(rm $names[$i]);

		#TODO: error handling if nucmer fails

		#Create the contigs array
		open(my $fhdelta, "<",$nucmer_tempy) 
		    or die "cannot open < ".$nucmer_tempy." : $!";
		my ($c1,$c2) = get_contigs($fhdelta);
		%contigs1 = %$c1;
		%contigs2 = %$c2;
		close($fhdelta);
		#qx(rm -f $nucmer_tempy."delta");


		#Check positions.
		if (not(($qi==0)and($ri==1))) {
		    if ($pre_comb[0]==$qi) {
			%contigs1 =  position (\%contigs_p1,\%contigs1);		    
		    }
		    elsif ($pre_comb[0]==$ri) {
			%contigs2 = position (\%contigs_p1,\%contigs2);;
		    }

		    if ($pre_comb[1]==$qi) {
			%contigs1 = position (\%contigs_p2,\%contigs1);;
		    }
		    elsif ($pre_comb[1]==$qi) {
			%contigs2 = position (\%contigs_p2,\%contigs2);;
		    }
		}
		
		#Writes the new subsets to the new temp files. Then changes the names in the filename array.
		$subQ_tempy = "sub_tempQ".$qi.$ri.".fasta"; 
		$subR_tempy = "sub_tempR".$qi.$ri.".fasta"; 
		filter($filenames[$qi],$subQ_tempy,\%contigs1);
		filter($filenames[$ri],$subR_tempy,\%contigs2);
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;

		#Sets the iformation needed in next iteration.
		@pre_comb = ($qi,$ri);
		%contigs_p1 = %contigs1;#Store information for position checks in next iteration
		%contigs_p2 = %contigs2;#Store information for position checks in next iteration
		$i++;
	    }
	}
    }

    #return tempname (@names);
    #qx(rm -f sub_temp*.fasta);
    my $test = makeFasta(\%contigs1,$subQ_tempy);

    return $test; 
}

sub tempname {
    for my $delta (@_) {
	#my ($a,$b) = get_contigs($delta);
    }

    return @_;
}

sub makeFasta {
    my %fast;
    my %contigs = %{$_[0]};
    my $filename = $_[1];
    my @keys = keys (%contigs);

    my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
    my $filtered = Bio::SeqIO->new(-file => ">finalpaoutput.fasta", -format => 'fasta' );
    while (my $seq = $fh->next_seq)
    {
	my $contigname = $seq->display_id;

	for my $cid (@keys)
	{
	    if ($cid eq $contigname)
	    {
		my @t = @{$contigs{$cid}};
		my $start = $t[0][0];
		my $end = $t[0][1];
		$fast{$cid} = [$start,$end];
		my $subseq = $seq->subseq($start,$end);
		$seq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname");
		$filtered->write_seq($seq);
	    }
	}
    }
    return \%fast;
}

#function that goes through the .delta file and stores the contig ids and their positions in an array.
sub get_contigs {
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

# Function that filters away contigs and creates subsets of original assembly files.
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

			    if ($cid eq $cidp) {
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
    }

    return %cpos = %cpos;
}

sub position2 {
    my %cpos;
    my %c1 = %{ $_[0] };
    my %c2 = %{ $_[1] };

    for my $cid1 (keys(%c1)) {
	for my $cid2 (keys(%c2)) {
	    if (defined($c1{$cid1})) {
		my @pos;
		for my $po_cur (@{ $c1{$cid1} }) {
		    if (defined($c2{$cid2})) {
			for my $po_pre (@{ $c2{$cid2} }) {

			    if($cid1 eq $cid2) {
				my $startP = min($po_pre->[0],$po_pre->[1]);
				my $endP = max($po_pre->[0],$po_pre->[1]);
				my $rangeP = Bio::Range->new(-start => $startP, -end => $endP);
				my $startC = min($po_cur->[0],$po_cur->[1]);
				my $endC = max($po_cur->[0],$po_cur->[1]);
				my $rangeC = Bio::Range->new(-start => $startC, -end => $endC);
				my ($start,$end,$strand) = $rangeC->intersection($rangeP);

				if(defined($start)) {
				    push (@pos ,[ $start,$end ]);
				    $cpos{$cid1} = [ @pos ];
				}
			    }
			}
			
		    }
		}
	    }
	}
    }

    return %cpos = %cpos;
}
1;
