package pairwiseAlignment;
use strict;
use warnings;
use lib do { # Incase bioperl is not installed using aptitude or apt-get. Library containing the bioperl library should be placed in subfolder relative to this file. The part "BioPerl-1.6.901/install/share/perl" bellow is that subfolders path relative to the pairwiseAlignment.pm file.
    use Cwd qw/realpath/;
    realpath($0) =~ /^(.*)\//;
    "$1/BioPerl-1.6.901/install/share/perl";
};
use Bio::Perl;
use Data::Dumper;
use Bio::SeqIO;

sub min { # Takes a any number of numbers and returns the smallest
    return (sort {$a<=>$b} @_)[0]; 
}

sub max { # Takes a any number of numbers and returns the biggest
    return (sort {$b<=>$a} @_)[0]; 
}

sub savey { # Helper function to save output to a textfile so you can check that stuff are correct. Not needed in finished program.
    my ($filename,$thingy) = @_;

    open(my $fhtest, ">",$filename) 
	or die "cannot open < ".$filename." : $!";
    print $fhtest "Thingy:\n";
    print $fhtest Dumper \$thingy;
    print $fhtest "\n";
    close($fhtest);
}

sub pa { # Main subroutine. Calls uses the other subroutines in pairwiseAlignment.pm to perform the iterative pa. It's final product is a fasta file called finalpaoutput.fasta
    my $fhdelta;
    #Temporary filenames for subset of assemblyfiles 
    my $subQ_tempy; 
    my $subR_tempy;
    #Temporary filename for deltafile
    my $nucmer_tempy;
    my (%contigs1,%contigs2); # Current contig hashes. Maps contig id to start and end position of matches in alignments
    my (%contigs_p1,%contigs_p2); # Previous contig hashes. Maps contig id to start and end position of matches in alignments
    #my @contigs_pos; # For storage of all contig hashes.
    my @filenames = @_;
    my @pre_comb;#keeps track of which assemblies were combined in the previous iteration
    my $i = 0; # Counter variable
    my @names; 
    my %r_ids; # Contig ids of removed contigs

    for my $qi (0..$#filenames) {
	for my $ri (0..$#filenames) {
	
	    if ($qi != $ri) {
		$nucmer_tempy = "unf_nucmer_out_temp".$qi."_".$ri; 
		$names[$i] = $nucmer_tempy.".delta";

		#TODO: error handling if nucmer fails
		qx(nucmer --prefix="$nucmer_tempy" "$filenames[$qi]" "$filenames[$ri]");
		$nucmer_tempy = "nucmer_out_temp".$qi."_".$ri.".delta"; 
		qx(delta-filter -i 99 -l 100 $names[$i] > "$nucmer_tempy");		
		qx(rm $names[$i]); # Cleanup of unfiltered deltafiles

		#Create the contigs array
		open(my $fhdelta, "<",$nucmer_tempy) 
		    or die "cannot open < ".$nucmer_tempy." : $!";
		my ($c1,$c2) = get_contigs($fhdelta);
		%contigs1 = %$c1;
		%contigs2 = %$c2;
		close($fhdelta);

		qx(rm -f $nucmer_tempy."delta"); # Done with the deltafiles. Cleanup.

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
		%r_ids = filter($filenames[$qi],$subQ_tempy,\%contigs1,\%r_ids,$qi,$i);
		%r_ids = filter($filenames[$ri],$subR_tempy,\%contigs2,\%r_ids,$ri,$i);
		$filenames[$qi] = $subQ_tempy;
		$filenames[$ri] = $subR_tempy;

		#Sets the iformation needed in next iteration.
		@pre_comb = ($qi,$ri);
		#$contigs_pos[$i] = \( \%contigs1, \%contigs2 );
		%contigs_p1 = %contigs1;#Store information for position checks in next iteration
		%contigs_p2 = %contigs2;#Store information for position checks in next iteration
		$i++;
	    }
	}
    }
    makeFasta(\%contigs1,$subQ_tempy);
    my $rem = '';
    for $a (keys(%r_ids)) {
	my @r = keys(%{$r_ids{$a}});
	$rem .= "\nNumber of removed contigs in assembly ".$a.": ".$#r."\n"; 
    }
    savey("removed_contigs.txt",\%r_ids);
    #qx(rm -f sub_temp*.fasta);# Cleanup

    return $rem; 
}


# Here the final result gets saved into a fasta file
sub makeFasta {

    my %contigs = %{$_[0]};
    my $filename = $_[1];
    my @keys = keys (%contigs);

    my $fh = Bio::SeqIO->new(-file => "$filename", -format => "fasta");
    my $filtered = Bio::SeqIO->new(-file => ">finalpaoutput3.genbank", -format => 'genbank' );
    while (my $seq = $fh->next_seq)
    {
	my $contigname = $seq->display_id;

	for my $cid (@keys)
	{


	    if ($cid eq $contigname)
	    {

		my @t = @{$contigs{$cid}};#Array of all the position hits of a given contig id

		    
		for my $pos (@t) {

		    my $start = $pos->[0];
		    my $end = $pos->[1];
		    my $subseq = $seq->subseq($start,$end);
		    my $seq2 = Bio::Seq->new(-seq => "$subseq", -id => "$contigname");
		    $filtered->write_seq($seq2);

		} 
	    }
	}
    }
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
	my ($filename, $output, $cont,$rem, $assembly,$align) = @_;
	my @contigs = keys (%{$cont});
	my %removed_ids = %{$rem};
	my %clist;
	if (defined($removed_ids{$assembly})) {
	    %clist = %{$removed_ids{$assembly}};
	}
	my $fh = Bio::SeqIO->new(-file => "$filename", 
				 -format => "fasta");
	my $filtered = Bio::SeqIO->new(-file => ">$output", 
				       -format => 'fasta' );
	while (my $seq = $fh->next_seq) {
		my $contigname = $seq->display_id;
		my $hits = 0;
		for my $cid (@contigs) {
			if ($cid eq $contigname) {

                                # This part is left out for now but can be useful later to see how it changes the result
				#my $subseq = $seq->subseq($contigs[$index][1],$contigs[$index][2]);
				#$subseq = Bio::Seq->new(-seq => "$subseq", -id => "$contigname"."\.$hits");
			
			    $hits++;
			    $filtered->write_seq($seq);
			}
		}
		if ($hits == 0) {		    
		    $clist{$contigname} = $align;
		    $removed_ids{$assembly} = \%clist;
		}
	}
	return %removed_ids;
}

sub position {
    my %cpos;
    my %c_p = %{ $_[0] };#contig information from previous iteration
    my %c   = %{ $_[1] };#contig information from current iteration
    my $minlength = 10;

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
				    if ((  $end-$start ) > $minlength ) {
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
    }

    return %cpos = %cpos;
}
1;
