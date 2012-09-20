package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
#use Exporter;


sub pa {
  for my $qi (0..($#_-1)) {
      for my $ri (($qi+1..$#_)) {
	  my $nucmer_tempy = "nucmer_out_temp".$qi."_".$ri; 
	  qx(nucmer --prefix="$nucmer_tempy" "$_[$qi]" "$_[$ri]");
          #get .delta file
          #delta-filter to remove alignments with low identity
          #.delta file -> fasta format
          #Extract consensus and other alignment information (length, identity)
          #Output extracted information to main/file/stdout
	  qx(rm -f "$nucmer_tempy.delta");
      }
  }
}
1



