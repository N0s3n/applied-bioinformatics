package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
#use Exporter;


sub pa {
    my $minid = 90;
    my $minlength = 20;
    my $alignment_p;
    my $alignment_c;

    for my $qi (0..$#_) {
	for my $ri (0..$#_) {
	    if ($qi not $ri) {

		if (($qi != 0) & ($ri != 1)) {
		    open($fh,'<',$_[$qi])
		    write($nucmer_tempy.delta);
		    #Filter out contigs
                    #get_info($nucmer_tempy);
		}

		my $nucmer_tempy = "nucmer_out_temp".$qi."_".$ri; 
		qx(nucmer --prefix="$nucmer_tempy" "$_[$qi]" "$_[$ri]");
		#A = AiAj
		#Ai = A_subi
		#Aj = A_subj
		#get .delta file
		#delta-filter to remove alignments with low identity
                #qx(delta-filter -i="$minid" -l="$minlength" "$nucmer_tempy.delta");
		#.delta file -> fasta format
		#Extract consensus and other alignment information (length, identity)
		#Output extracted information to main/file/stdout
		qx(rm -f "$nucmer_tempy.delta");
	    }
	}
    }
}

sub get_info {
    #get position and contig information from .delta file using regular expressions
}
1



