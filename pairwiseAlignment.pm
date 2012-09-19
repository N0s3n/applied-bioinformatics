package pairwiseAlignment;
#use base qw( Exporter );
#our @EXPORT_OK = qw(pa);
use strict;
use warnings;
use Exporter;

sub pa {
  my $variable;
  for my $qi (0..($#_-1)) {
      for my $ri (($qi+1..$#_)) {
	  #$variable .= $_[$qi].$_[$ri];
	  qx(nucmer $_[$qi] $_[$ri]);
      }
  }
  return $variable;
}
1



