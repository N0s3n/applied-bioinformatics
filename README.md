<h3>name</h3>
<p>CSF - Common Sequence Finder</p>

<h3>synopsis</h3> 

<p>perl main.pl [FILE] </p3>

<h3>description</h3>

<p>CSF is a tool written in perl as a student project. The tool will take assemblies created from the same sequencing data and do a multiple set of tasks.
</p>
<p>
The first thing it'll do is to calculate some basic statistics on the assemblies themselves, such as: Total length, Number Of contigs, N50, N90 and the Largest contig.
</p>
<p>
The tool will try to come up with a subset of contigs that are represented in all the the assemblies. It has been implemented in two ways. 
The first is to run the assemblies in a multiple alignment program called Mauve and extract the common contigs from the output. 
The latter is do to pairwise alignments between two assemblies and take the intersection between them, if there are more than two assemblies it will rerun the pairwise alignment, between the intersection and the next assembly. Repeats until all the assemblies have been processed.
There arose a fork on the pairwise alignment implementation. Both of them are present in the tool and both produce an output. 
</p>
<p>
In theory all three should produce the same output. But in practice that is not always the case.
</p>
<p>
An attempt to extract difficult regions is also implemented. Mauve produces a file that contains coordinates where different contigs from different assemblies matched. 
If a contig only has coordinates represented in one assembly then it is a difficult region. 
This has not been fully implemented. The tool produces sequences that are present in only one assembly. But Mauve appends several contigs and sets the coordinates from where the first contig starts and the last contig ends. The thing that is not implemented is to say from what contig/contigs the sequence came from. Currently it only says from which Assembly the difficult region is from. 
</p>
<h3>requirements</h3>
<p>Perl</p>

<p>
Perl-modules that have to be installed:<br>
File::Slurp<br>
Text::Table<br>
Bio::Perl<br>
</p>
<h3>files</h3>
<p>
Input:
Currently the only files that are necessary are the assembly files in fasta format.
ProgressiveMauve is picky when reading the fasta files. The fasta header can not include white spaces, and the sequences can't be in only one line. 
</p>
<p>
Output:<br>
Saves files to the output directory:
statistics.txt
</p>
<p>
Mauve files:<br>
commonContigs.fasta  <br>             
uniqueContigs.fasta
</p>
<p>
pairwiseAlignment:<br>
finalpaout.fasta<br>
finalpaout_removed_contigs.txt<br>
finalpaout_removed_contigs_total.txt<br>
</p>
<p>
pairwisebioperl:<br>
finalpairwisebioperl.fasta
</p>

<h3>diagnostics</h3>
<p>
If one of the assembly files have the improper structure stats.pm(statistics) will produce an error message. 
</p>
<h3>author</h3>
Written by Susanna Trollvad, Björn Viklund, Anders Sjölander and Anders Gustafsson.

--------------------------------------------------------------------------------
<p>
subroutine overview:

stats(statistics) (String -> filename) -> (String -> "Assembly \t Length \t Number of contigs \t N50 \t N90 \t Largest Contig") <br>

mauveCleaner(mauveParser) (String -> Alignment from mauve including information about the alignment) -> (String -> Multiple alignment)<br>
mauveCleaner(getAlignments) (String -> Multiple alignment, Scalar -> Number of assemblies) -> (Array -> Contig sequences that are represented in all assemblies)<br>

backbone(unmatchedCoo) () -> (AoAoA -> Coordinates from the .backbone file where only one contig is present)<br>
backbone(getPos) (String -> filename) -> (Array -> Starting positions of each contig)<br>
backbone(collapseAssembly) (String- > filename) -> (String -> All contigs on one line with no headers)<br>
backbone(getSubString) (String -> Collapsed sequence, Array -> Coordinates)

---

pairwiseAlignment(pa) (Array reference -> filenames,Hash reference -> options_pa(optional parameter)) -> (String -> "final file name")<br>
pairwiseAlignment(min) (Array -> integers) -> (Int -> smallest number)<br>
pairwiseAlignment(min) (Array -> integers) -> (Int -> largest number)<br>
pairwiseAlignment(savey) (String -> filename,Scalar -> content) -> (void)<br>
pairwiseAlignment(makeFinalFile) (Hash reference -> contigids,String filenamein,String -> finalfilename,String -> format) -> (String -> filenameout)<br>
pairwiseAlignment(get_contigs) (Hash reference -> contigs1,Hash reference -> contigs2) -> (Hash reference -> contigs1,Hash reference contigs2)<br>
pairwiseAlignment(filter) (String -> filename,String -> filenameout,Hash reference -> ) -> (Int -> smallest number)<br>
pairwiseAlignment(position) (Hash reference -> contigs previous, Hash reference -> contigs current)-> (Hash -> contigs)<br>

---

pairwisebioperl(get_contigs) (String -> nucmer output name) -> (Hash -> contigid keys start and end position array values)<br>
pairwisebioperl(filter) (Hash -> contigid keys start and end position array values, String -> name of output, String -> name of input fasta files) -> ()<br>


--------------------------------------------------------------------------------

Future work:<br>

Mauve part:<br>
Currently when parsing the mauve output and obtaining the sequences, only sequences that exists in all assemblies are saved. A good thing to implement here is to be able to choose if the contig have to exists in all or maybe all except one assembly.<br>

When parsing the .backbone file. Mauve concatenates contigs that don't match and returns the coordinates from the first contigs first position and the last contigs last position. A thing to implement here is a function that returns from which assembly and which contigs those coordinates came from.<br>

A small script that checks the assembly fasta files for improper headers, and also check that the sequences are not written in one line. For instance about 60-80 nucleotides per line.<br>

Right now progressiveMauve runs its default values. A good thing to implement is some parameter that lets the user choose if he or she want some specific parameters that progressiveMauve should run.<br>

The .backbone file can be used to extract coordinates where every contig is present. Would maybe be a fun thing to check if the output from the multiple alignment and the coordinates from the .backbone file is equal.<br>

Clean the commonContigs.fasta file from gaps.  <br><br>

pairwisebioperl:<br>
Perhaps the most obvious problem is the simplicity of the approach. As is, the subroutine simply discards every bit of sequence data not found in all assemblies. The most obvious alternative would be to store the number of assemblies the sequence is found in.<br>

The way that the contigs are progressively renamed can be somewhat problematic. A different way to distinguish between subsequences originating from the same contig may be desirable. Simply continually updating the start and end position data rather than creating entirely new fasta files would allow this, as well as rendering backtracking far easier.<br>

There are still a few mismatches when running assemblies against themselves. Further experimentation with thresholds may improve the results.<br>

The subroutine has not been made with scaffolds in mind, does not deal with them well, and further work is needed to handle them and their gaps.<br><br>

pairwiseAligment:<br>

The final output fasta file contains the sequences common to all the assemblies. But the ids of these sequences are the ids of the contigs in one of the last temporary subset fasta files. It should be changed so that when you create the final fasta file the sequence is assigned a new unique identifier. <br>

It would be very useful if the algorithm kept track of which contigs in the assemblies contributed to the found common sequence in the end. One way of doing this could be to write it directly into the temporary subset fasta files in the filter function. However as it is now filter does not receive enough information to do this fully. The contig hashes for the query and reference sequences does not contain information connecting the contigs which is needed in order to keep track of the contributing contigs.<br>

The module picks out the sequences common to all assemblies, but a good extension of this could be to make it so that it can also pick out sequences that are represented in most assemblies. To do this you could for example use a threshold value and after each alignment make a check if a certain contig has failed to make a match or make an intersecting match above this threshold, and then filter it away first after this.<br>

In the position subroutine of pairwiseAlignment you check if there is an intersection between matches on contig A in the current iteration and the matches against contig A in previous iteration. This makes the algorithm very strict and will remove a lot of the contigs from the assembly. A way to change this would be to remove the position function calls from the pa subroutine altogether.
</p>
