#package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
use Exporter;
#use File::Slurp;

#my %paired = pa("out.delta");
#my @files = ("GMBLNrPeC02_ctg.fasta","GMBLNrPeM01S_ctg.fasta");
#my $postfilter;
#my $out;
#my @outfiles;
#foreach (@files)
#{
#	open ($out, ">", "$_" . ".output");
#	$postfilter .= filter($_, %paired);
#	print $out $postfilter;
#	close $out;
#	$postfilter = "";
#}
sub pa {
    my $minid = 90;
    my $minlength = 20;
    my $alignment_p;
    my $alignment_c;
    my $fh1;
    my $fh2;
    my $nucmer_tempy;
    my %filtered;
    my $contigs;


    for my $qi (0..$#_) {
	for my $ri (0..$#_) {

	    #OY! Temporary for testing
	    open(my $fh1, "<",$_[$qi]) 
		or die "cannot open < ".$_[$qi].": $!";
		    #open ($fh2,'<',$nucmer_tempy.".delta");
		    #$assembly = read_file($_[$qi]);
		    #$delta = read_file($nucmer_tempy.delta);

	    %filtered = get_contigs($fh1);
	
	    #Filter out contigs
	    #get_info($nucmer_tempy);
	    if ($qi != $ri) {

		if (($qi != 0) & ($ri != 1)) {
		    open(my $fh1, "<",$_[$qi]) 
			or die "cannot open < ".$_[$qi].": $!";
		    #open(my $fh2, "<",$nucmer_tempy.".delta") 
#			or die "cannot open < ".$nucmer_tempy.".delta".": $!";
		    
		    #$contigs = get_contigs($fh1);
	
		    #Filter out contigs
                    #get_info($nucmer_tempy);
		}

		#my $nucmer_tempy = "nucmer_out_temp".$qi."_".$ri; 
		#qx(nucmer --prefix="$nucmer_tempy" "$_[$qi]" "$_[$ri]");
		#A = AiAj
		#Ai = A_subi
		#Aj = A_subj
		#get .delta file
		#delta-filter to remove alignments with low identity
                #qx(delta-filter -i="$minid" -l="$minlength" "$nucmer_tempy.delta");
		#.delta file -> fasta format
		#Extract consensus and other alignment information (length, identity)
		#Output extracted information to main/file/stdout
	    }
	}
    }
    #qx(rm -f "$nucmer_tempy.delta");
    #qx(rm -f "$nucmer_tempy.txt");
    return %filtered;
}

#function that goes through the .delta file and stores the contig ids and their positions in a hash.
sub get_contigs {
    my %contigs;
    my $c1;
    my $c2;

    for my $fh (@_) {
	for my $line (<$fh>) {
	    if ($line =~ /^>([^\s]+)\s+([^\s]+)/) {
		$c1 = $1;
		$c2 = $2;
	    }
	    if ($line =~ /^(\d+)\s(\d+)\s(\d+)\s(\d+)\s\d+\s\d+\s\d+$/) {
		$contigs{$c1} = ($1,$2);
		$contigs{$c2} = ($3,$4);
	    }
	}
    }
    return %contigs;
}

sub filter {
	my ($filename, %contigs) = @_;
	my @keys = keys %contigs;
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
				}
			}
		}
		elsif ($found == 1)
		{
			$content .= $in;
		}
	}
	return $content;
}
1