#!/usr/bin/perl
use warnings;
use strict;
use Bio::SimpleAlign;
use Bio::AlignIO;

my $in = Bio::AlignIO->new(-file => "2.txt", -format => "fasta");
open (my $out, ">", "aligned2.txt") or die "Do you want your possesions identified?";
while (my $aln = $in->next_aln)
{
	my $iupac_consensus = $aln->consensus_iupac(); #Creates a consensus sequence. Gaps in any of the sequences result in non-capital letters.
	#$iupac_consensus =~ tr/[a-z]/-/; #replaces gaps with -
	#$iupac_consensus =~ s/[a-z]//g; #rips out gaps
	#$iupac_consensus =~ tr/[a-z]/[A-Z]/; #hides the fact that there ever were gaps. Utilize to confuse people.
	#print $out $iupac_consensus;
	print $iupac_consensus;
}
