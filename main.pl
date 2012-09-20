use strict;
use warnings;
use Getopt::Long;
use BIO::AlignIO;

use mauveCleaner qw ( mauveParser );
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
#	print $_;
}

#Run mauve
my $concatFilenames = join (" ", @filenames);
my $mauve = `progressiveMauve $concatFilenames`;
$mauve = mauveParser($mauve); 

print "$mauve";

# PA
