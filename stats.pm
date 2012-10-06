use strict;
use warnings;
use Exporter;

#The subroutine statistics will get a filename as input and will return the:
#"Assembly", "Length", "NuOfCon","N50", "N90", "Largest Contig"
sub statistics {
  my ($filename) = @_;
  open my $fh,$filename or die "$!,$?";
  my @array;
  my $count=-1;
  my $totLength = 0;
  my $numOfContigs = 0;
  while (my $line = <$fh>) {
    if ($line =~ /^>/) {
      $count++
    }
    else {
      $array[$count] += length($line)-1; 
      $totLength += length($line)-1;
    }
  }

#To get the N50 and N90 the @array that contains all the sequence lengths are sorted from the largest to the smallest.
#If the sum of the contig lengths > total length * 0.5, then the current contig length is stored in a scalar.
#And if the sum of the contig lengths > total length * 0.9, then it has found all the values and will return the statistics.
	my @sort = sort{$b <=> $a} @array;
	my $largestContig = $sort[0];
	my $n50;
  my $n50ret = 0;
	foreach my $val (@sort) {
		$n50 += $val;
		if ($n50 >=  $totLength*50*0.01 && $n50ret == 0) {
			  $n50ret = $val; 
		}
    elsif($n50 >=  $totLength*90*0.01) {
      chomp $filename;
      $count++;
      my $ret = "$filename \t $totLength \t $count\t  $n50ret \t $val\t $largestContig \n";
      return $ret;
    }
	}
}
1
