use strict;
use warnings;
use Exporter;


sub mauveParser {
  my ($mauveInput) = @_;
  open my $fh,"<", \$mauveInput or die $!;

  my @alignment;
  my $bool = 0;
  while (<$fh>) {
    if (/^>/) {
      $bool = 1;
    }
    if ($bool == 1) {
      push(@alignment,$_);
    }
  }
  my $finalAlignment = join("",@alignment);
  return $finalAlignment;
}
1
