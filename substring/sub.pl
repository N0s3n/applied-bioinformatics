use strict;
use warnings;

open my $fh, "<", "contig.fasta" or die $!;
my @positionArray;
my $coordinate = 1;

while (<$fh>) { 
  if (/^>/) {
    push(@positionArray, $coordinate);
  }
  else {
    $coordinate += (length($_)-1);
  }
}
close $fh;

open $fh, "<", "contig.fasta" or die $!;
my $joinedString = "";
while (<$fh>) {
  if ($_ !~/^>/) {
    chomp;
    $joinedString = join("",$joinedString, $_);  
  }
}
close $fh;

my @cordArray;
open my $backbonefh, "<", ".backbone" or die "pallar inte livet, $!";
while (<$backbonefh>) {
  my @line = split(/\t/);
  push(@cordArray, [$line[0], $line[1]]);
}
#$cordArray[1][1]" 


my @substring;
foreach (@cordArray) {
  push(@substring, substr($joinedString,$_->[0]-1,1+($_->[1]-$_->[0])));
}
foreach (@substring) { 
  print $_ . "\n";

}


