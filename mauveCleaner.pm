use strict;
use warnings;
use Exporter;

#mauveParser cleans the output from mauve. The first ~100 lines is only mauve saying how much of the alignment it has done. It returns only the multi fasta alignment. 
sub mauveParser {
	my ($mauveInput) = @_;
	open my $fh,"<", \$mauveInput or die $!;

	my @alignment;
	my $bool = 0;
	while (<$fh>) {
	  if (/^>/) {
		$bool = 1;
	  }
	  if ($bool == 1) {
		push(@alignment,$_);
	  }
	}
	my $finalAlignment = join("",@alignment);
	return $finalAlignment;
}


#getAlignments will loop through the output from mauve and extract only the contigs that are represented in all the assemblies. It will skip alignments that contains 0-0 coordinates. And if not all the assemblies are present.
#To check if not all the assemblies are present an if statement will check if the #loopCount is equal to the assembly number contained in the regexp. If one assembly is missing those numbers won't match and won't be pushed. 
#To skip an alignment the script will zero the $headercount variable and it will never reach the if statement where it actually get pushed into the @alignmentArray. 
sub getAlignments {

	my $headerCount = 0;
	my @alignmentArray;
	my $fasta = ">";
	my $loopCount = 1;
	my $sequence = "";
	my ($alignmentString, $noOfHeaders) = @_;
	open my $fh,"<", \$alignmentString or die $!;
	  
	while (<$fh>) {
		#print "$loopCount \n";
		if (/^>\s(\d+):(\d+)-(\d+)/) {
			$fasta .= "$2-$3 ";
			if (/0-0/) {
			  $headerCount = 0;
			}
			elsif ($1 != $loopCount) {
			  $headerCount = 0;
			}
			$loopCount++;
			$headerCount++;
			$sequence = "";
		}

		elsif (/^=/) {
			if ($headerCount == $noOfHeaders) {
				$fasta .="\n$sequence";
				push (@alignmentArray,$fasta);
		  	}
			$headerCount = 0;
			$fasta = ">";
			$loopCount = 1;
		}

			else {
			  #$sequence =~ s/\-//gs,$sequence;
			  $sequence .= $_;
			  #print $sequence;
			  #print  "\n";
			  #$sequence .= s/-//g;
			}
	}
	return @alignmentArray;
}
1
