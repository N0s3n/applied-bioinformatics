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
  push(@substring, substr($joinedString,abs($_->[0])-1,1+abs($_->[1]-$_->[0])));
}
foreach (@substring) { 
#  print $_ . "\n";

}


#parse .backbone

open  $fh , "../.backbone" or die "$!";
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
    #  push(@outArray[$indices[0]/2], ($xs[$indices[0]],$xs[$indices[1]]));
    @array = [$indices[0], $indices[1]];
    push (@outArray -> ($indices[0]/2) , @array);
    #push (
    #print $outArray[0][0][0] . "\n";
    #print $array[0]->[0];
  }
  print $outArray[0][1][0] . "\n";



  #print "$nzeroes \n";
  #my $print = join(' ', @indices);
  #print $print . "\n";
}  
