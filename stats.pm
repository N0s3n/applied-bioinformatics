use strict;
use warnings;
use Exporter;
sub statistics {
	my ($filename,$nvalue) = @_;
	open my $fh,$filename or die "$!,$?";
	my @array;
	my $count=-1;
	my $totLength = 0;
	while (my $line = <$fh>) {
		if ($line =~ /^>/) {
			$count++
		}
		else {
			$array[$count] += length($line)-1; 
		 	$totLength += length($line)-1;
		}
	}

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
      my $ret = "$filename \t $totLength \t  $n50ret \t $val\t $largestContig \n";
      return $ret;
    }
	}
}
1
