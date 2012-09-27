use warnings;
use strict;
use filter;

my %paired = pa("out.delta");
my @files = ("GMBLNrPeC02_ctg.fasta","GMBLNrPeM01S_ctg.fasta");
my $postfilter;
my $out;
foreach (@files)
{
	open ($out, ">", "$_" . ".output");
	$postfilter .= filter($_, %paired);
	print $out $postfilter;
	close $out;
	$postfilter = "";
}
