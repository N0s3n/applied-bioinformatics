use strict;
use warnings;

my @filenames = @ARGV;
#PA
#Send all combinations of filenames to NUCmer
my $i = 0;
my $j = 1;
my $arrlength = $#ARGV;
while ($i < $arrlength){
	while ($j <= $arrlength){
		#print "$filenames[$i] $filenames[$j] \n";
		system("nucmer -p $i$j $filenames[$i] $filenames[$j]"); 
		# the -p designates output filename, could probably be polished.
		# any extra options for nucmer? output files don't look great as it is.
		$j++;
	}
	$i++;
	$j = $i + 1;
}
