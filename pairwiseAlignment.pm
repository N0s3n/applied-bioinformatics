package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
#use Exporter;
#use File::Slurp;

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

sub get_contigs {
    #my @contigs;
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
		$contigs{$c1} = " c1: ".$1."\n";
		$contigs{$c2} = " c2: ".$3."\n";
	    }
	}
    }
    return %contigs;
}
1
