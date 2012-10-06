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
#my $nvalue = 50;
#GetOptions ("nvalue=i" => \$nvalue);
my @filenames = @ARGV;

#Create output folder
mkdir "output";

#Basic statistics
my @statsArray;
foreach (@filenames) {	
	my $statsOutput = statistics($_);
	push (@statsArray,$statsOutput);
}

#Print statistics to output/statistics.txt
my $statTable = Text::Table->new("Assembly", "Length", "NuOfCon","N50", "N90", "Largest Contig");
$statTable->load(@statsArray);
#open my $myfile, ">output/statistics.txt";
#print $myfile $statTable;
#close $myfile;
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

#Print common contigs to file
my $commonOutput = join "", @commonContigs; 
#Clean spaces and gaps from sequence
#my $test = '';
#open my $fh,"<", \$commonOutput or die $!;
#while (<$fh>) {
#  if (/^>/) {
#     $test.=$_;
#  }
#  else
#}
open my $myfile, ">output/commonContigs.fasta";
print $myfile $commonOutput;
close $myfile;
#print $commonOutput;


#Parse the .backbone file
#--------------------------------------------------------------------------------

#Get positions of each file
#Define variables
my @assemblyCtgCoo;
my @chompedSeq;
my $currentAssembly = 0;
my @difficultContigs;
my @position;

#Store the coordinates of each unmatched sequence. The array is an array in an array in an array. 
#The first level is which assembly the coordinate comes from. The second level is which unmatched sequence. And the third level is either the start or the end of the coordinate. 
my @unmatchedContigs = unmatchedCoo(); 

foreach my $filename (@filenames) {
    #getPos will store all the starting positions of the contigs in an assembly.
    my @coordinates = getPos($filename);
    push @assemblyCtgCoo, \@coordinates;

    #This was an attempt to map specific unmatched contigs back to which assembly and position it came from. But mauve seemed to append over several contigs. When we figured that out, the course was near its end and we did not have time to implement it properly.
    foreach my $unmatched (@{$unmatchedContigs[$currentAssembly]}) {
      my @tempArray = grep @{$assemblyCtgCoo[$currentAssembly]}[$_] >= $unmatched->[0], 0 .. $#{$assemblyCtgCoo[$currentAssembly]};
      push @{$position[$currentAssembly]}, $tempArray[0]; 
    }
    
    #collapseAssembly chomps each assembly so that we're working on a single long sequence to be able to extract the sequences from the coordinates recieved from the .backbone file. getSubString will get the chomped sequence and all the coordinates and return the different difficult "regions" in an array. 
    push @chompedSeq, collapseAssembly($filename);
    @{$difficultContigs[$currentAssembly]} =  getSubString($chompedSeq[$currentAssembly], @{$unmatchedContigs[$currentAssembly]});
    $currentAssembly++;
}

#Extracts the unique contigs and stores it in fasta format in one single variable that is later written to output/uniqueContigs.fasta
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

#Run statistics on the created files.
$statTable->load(statistics("output/commonContigs.fasta"));
$statTable->load(statistics("output/uniqueContigs.fasta"));
open  $myfile, ">output/statistics.txt";
print $myfile $statTable;
close $myfile;


# PA
my %pa_options;# This is an optional second parameter to pa(). If given it must contain:
# keys: ('fn','minl','mi','format') (string,int,float,string)
# fn = filename, with no format identifier, eg. 'finalfile'
# minl = minimum length which will take away intersection(found in subroutine position() in pairwiseAlignment.pm) two matched positions in the same contig
# mi, minimum % identity between matched contigs 
my $removed_ctg = pa(\@filenames);
