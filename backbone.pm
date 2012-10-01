#!/usr/bin/perl
use warnings;
use strict;
use Exporter;


# Will return the starting positions of each contig in a file.
sub getPos {
	my ($filename) = @_;
	open my $fh, "<", $filename or die $!; 
	my @positionArray;
	my $coordinate = 1;

	#For each fasta header, the position will be pushed into @positionArray 
	while (<$fh>) { 
	  if (/^>/) {
		push(@positionArray, $coordinate);
	  }
	  else {
		$coordinate += (length($_)-1);
	  }
	}
	close $fh;
	return @positionArray;
}

#Will return whole assembly on one line.
sub	collapseAssembly {
	my ($filename) = @_;
	open my $fh, "<", $filename or die $!;
	my $joinedString = "";
	#Every single line that is not a header will be joined.
	while (<$fh>) {
	  if ($_ !~/^>/) {
		chomp;
		#$joinedString = join("",$joinedString, $_);
		$joinedString .= $_;  
	  }
	}
	close $fh;
	return $joinedString;
}

#Parse the .backbone file that mauve returns. Will return the coordinates of each unmatched sequence in all assemblies. And keeps track of which coordinates goes to which assembly.  
sub	unmatchedCoo {
	open my $fh , ".backbone" or die "$!";
	my $nzeroes;
	my @indices;
	my @outArray;
	my @array;
	while (<$fh>) {
	  next if /^seq0/;
	  my @xs = split(/\s/);
	  $nzeroes = grep $_ == 0, @xs;
	  @indices = grep $xs[$_]!= 0, 0 .. $#xs;
	  if (scalar(@indices) == 2 ) {
		my @tusn = ($indices[0],$indices[1]);
		push @{$outArray[$indices[0]/2]}, [$xs[$indices[0]],$xs[$indices[1]]];
	  }
	}
	return @outArray;	
}

#Given the collapsed contigs from one assembly and the coordinates from a unique contig. The contig will be saved in @substrings. The coordinates will be mapped back to the sequence and the unmatched sequence will be returned.
sub getSubString {
	my ($joinedString, @coordArray) = @_;
	my @substrings;
	#print $joinedString;
	foreach my $cooPair (@coordArray) {	
		push(@substrings, substr($joinedString,abs($cooPair->[0])-1,1+abs($cooPair->[1]-$cooPair->[0])));
	}
	
	return @substrings;
	#return [1,0];
}



1;
