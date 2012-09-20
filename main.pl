use strict;
use warnings;
use Getopt::Long;
use BIO::AlignIO;

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

# MSA
#     my $multipleAlignments = qx(progressiveMauve inpara);
# loop
# fastaExtract.pl
# consensus module
# end loop

# PA
# Send all combinations of filenames to NUCmer
#my $variable = pa();
#print $variable . "\n";
