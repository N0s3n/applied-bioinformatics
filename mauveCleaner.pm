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



sub getAlignments {

  my $headerCount = 0;
  my @alignmentArray;
  my $fasta = ">";
  my $loopCount = 1;
  my $sequence = "";
  my ($alignmentString, $noOfHeaders) = @_;
  open my $fh,"<", \$alignmentString or die $!;
    
  while (<$fh>) {
    #print "$loopCount \n";
    if (/^>\s(\d+):(\d+)-(\d+)/) {
      $fasta .= "$2-$3 ";
      if (/0-0/) {
        $headerCount = 0;
      }
      elsif ($1 != $loopCount) {
        $headerCount = 0;
      }
    $loopCount++;
    $headerCount++;
    $sequence = "";
    }

    elsif (/^=/) {
      if ($headerCount == $noOfHeaders) {
        $fasta .="\n$sequence";
        push (@alignmentArray,$fasta);
      }
      $headerCount = 0;
      $fasta = ">";
      $loopCount = 1;
    }

    else {
      $sequence .= $_;
    }
  }
  return @alignmentArray;
}
#sub fgetAlignments {
#  my ($alignmentString, $noOfHeaders) = @_;
#  my @alignmentArray;
#  my $fasta = "";
#  my $headerCount = 0;
#  open my $fh,"<", \$alignmentString or die $!;
#  while (<$fh>) {
#    if (/=/) {
#      if ($headerCount == $noOfHeaders) {
#        push (@alignmentArray,$fasta);
#      }
#      $headerCount = 0;
#      $fasta = "";
#      next;
#    }
#    else {
#      if ($_ =~ /^>/) {
#        $headerCount++;
#      }
#    }
#    $fasta = join("",$fasta,$_);
#    $fasta .= $_;
#  }
#  return @alignmentArray;
#}
1
