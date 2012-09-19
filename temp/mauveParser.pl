use strict;
use warnings;



open my $fh, $ARGV[0] or die $!;

my $alignment = "";
my $start = 0;
while (<$fh>) {
  if ($_ =~ /^>/) {
    $start = 1;
  }
  if ($start == 1) {
    $alignment = join("",$alignment,$_);
  }
  print $start;
}
print $alignment;

