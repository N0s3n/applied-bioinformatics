package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
use Exporter;

sub pa {
  for my $qi (0..($#_-1)) {
      for my $ri (($qi+1..$#_)) {
	  qx(nucmer --prefix=temp $_[$qi] $_[$ri]);
      }
  }
}
1



