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
	my $bla = $totLength*$nvalue*0.01;
	foreach my $val (@sort) {
		$n50 += $val;
		if ($n50 >=  $totLength*$nvalue*0.01) {
			chomp $filename;
			my $ret = "$filename \t Total length: $totLength \t n$nvalue:$val \t largest contig: $largestContig \t Sum: $n50 \tMåste överstiga: $bla\n";
			return $ret;
			last;
		}
	}
	
}
1
