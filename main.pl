use strict;
use warnings;
use Getopt::Long;
use File::Slurp;
use Data::Dumper;
use Text::Table;
#use BIO::AlignIO;


#Import perl-modules
use mauveCleaner qw ( mauveParser getAlignments);
use pairwiseAlignment qw( pa );
use stats qw ( statistics );
use consensus qw ( mkConsensus );
use backbone qw ( getPos collapseAssembly unmatchedCoo );

#Declare variables
my $nvalue = 50;
my @filenames = @ARGV;
GetOptions ("nvalue=i" => \$nvalue);

#Create output folder
mkdir "output";

#Basic statistics
my @statsArray;
foreach (@filenames) {	
	my $statsOutput = statistics($_,$nvalue);
	push (@statsArray,$statsOutput);
}

#Print statistics
my $statTable = Text::Table->new("Assembly", "Length", "NuOfCon","N50", "N90", "Largest Contig");
$statTable->load(@statsArray);
open my $myfile, ">output/statistics.txt";
print $myfile $statTable;
close $myfile;
#print $statTable;

#Run mauve
#my $concatFilenames = join (" ", @filenames);
#my $mauve = `progressiveMauve $concatFilenames`;
#$mauve = mauveParser($mauve); 
#print "$mauve";

#Extract the common contigs from mauve output
my $fil = read_file("output.txt");
my @commonContigs = getAlignments($fil,scalar(@filenames));
#my @commonContigs = getAlignments($mauve,scalar(@filenames));

#Print contigs to file
my $commonOutput = join "", @commonContigs; 
open $myfile, ">output/commonContigs.fasta";
print $myfile $commonOutput;
close $myfile;
#print $commonOutput;


#Run statistics on the commonContigs
print statistics("output/commonContigs.fasta");
$statTable->load(statistics("output/commonContigs.fasta"));
open my $myfile, ">>output/statistics.txt";
print $myfile $statTable;
close $myfile;

#Parse the .backbone file
#--------------------------------------------------------------------------------

#Get positions of each file
my @assemblyCtgCoo;
my @chompedSeq;
my $currentAssembly = 0;
my @difficultContigs;
my @unmatchedContigs = unmatchedCoo(); 
my @position;


#print Dumper \$unmatchedContigs[0];
#print "Unmatched contigs:\n";
for (my $i = 0; $i <3; $i++) {
  for (my $y = 0; $y < scalar(@{$unmatchedContigs[$i]}); $y++) {
    #print join ("\t", @{$unmatchedContigs[$i][$y]});
    #print "\t\t";
  }
  #print "\n\n";
}
#print "Coordinates:\n";


foreach my $filename (@filenames) {
    my @coordinates = getPos($filename);
    push @assemblyCtgCoo, \@coordinates;
 
    #print join ("\t", @{$assemblyCtgCoo[$currentAssembly]});
    #print "\n";
    #Detta är sjölanders ansvar.
      foreach my $unmatched (@{$unmatchedContigs[$currentAssembly]}) {
        my @tempArray = grep @{$assemblyCtgCoo[$currentAssembly]}[$_] >= $unmatched->[0], 0 .. $#{$assemblyCtgCoo[$currentAssembly]};
        push @{$position[$currentAssembly]}, $tempArray[0]; 
        #last
      #print "$#{$assemblyCtgCoo[$currentAssembly]} \t $unmatched->[0]\n";
    }
    

    #print join ("\t", @position);
    #print "\n";
    push @chompedSeq, collapseAssembly($filename);
    @{$difficultContigs[$currentAssembly]} =  getSubString($chompedSeq[$currentAssembly], @{$unmatchedContigs[$currentAssembly]});
    $currentAssembly++;
}
my $firstcount = 0;
my $secondcount = 0;
my $uniqueContigs = '';
foreach my $firstline (@difficultContigs) {
  foreach (@{$firstline}) {
    $uniqueContigs .=  ">Assembly:$firstcount\t Contig:$secondcount \n$_\n";
    $secondcount++;
  }
  $secondcount = 0;
  $firstcount++;
}
open $myfile, ">output/uniqueContigs.fasta";
print $myfile $uniqueContigs;
close $myfile;



#print "\nPositions:\n";
#foreach my $pos (@position) {
#  for (my $i = 0; $i < scalar(@{$pos}); $i++) {
#      #print $pos->[$i]; 
#      #print "\t";
#  }
  #print "\n";
#}
#print Dumper @position;

#print $difficultContigs[1][1];

# PA
my %pa_options;# This is an optional second parameter to pa(). If given it must contain:
# keys: ('fn','minl','mi','format') (string,int,float,string)
# fn = filename, with no format identifier, eg. 'finalfile'
# minl = minimum length which will take away intersection(found in subroutine position() in pairwiseAlignment.pm) two matched positions in the same contig
# mi, minimum % identity between matched contigs 
my $removed_ctg = pa(\@filenames);
