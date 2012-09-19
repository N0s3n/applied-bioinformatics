use strict;
use warnings;
use Getopt::Long;

use pairwiseAlignment qw( pa );
use stats qw ( statistics );

my $nvalue = 50;
my @filenames = @ARGV;
GetOptions ("nvalue=i" => \$nvalue);

#Basic statistics
my @statsArray;
foreach (@filenames) {	
	my $statsOutput = statistics($_,$nvalue);
	push (@statsArray,$statsOutput);
}

#Print statistics
foreach (@statsArray) {
	print $_;
}


