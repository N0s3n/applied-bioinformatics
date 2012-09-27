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


#getAlignments need alot of optimization. Now it takes 1 min to run.
sub getAlignments {
  my ($alignmentString, $noOfHeaders) = @_;
  my @alignmentArray;
  my $fasta = "";
  my $headerCount = 0;
  open my $fh,"<", \$alignmentString or die $!;
  while (<$fh>) {
    if (/=/) {
      if ($headerCount == $noOfHeaders) {
        push (@alignmentArray,$fasta);
      }
      $headerCount = 0;
      $fasta = "";
      next;
    }
    else {
      if ($_ =~ /^>/) {
        $headerCount++;
      }
    }
    $fasta = join("",$fasta,$_);
  }
  return @alignmentArray;
}
1
