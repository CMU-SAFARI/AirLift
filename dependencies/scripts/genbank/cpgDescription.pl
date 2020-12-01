#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $argc = scalar(@ARGV);

if ($argc != 1) {
  printf STDERR "usage: cpgDescription.pl <buildDirectory>\n";
  printf STDERR "where <buildDirectory> is where the assembly hub is under construction.\n";
  exit 255;
}

my $wrkDir = shift;

if ( ! -d $wrkDir ) {
  printf STDERR "ERROR: can not find work directory:\n\t'%s'\n", $wrkDir;
  exit 255;
}

chdir $wrkDir;
# print STDERR `ls`;
my $accAsm = basename($wrkDir);
exit 0 if ( ! ( -s "bbi/${accAsm}.cpgIslandExtUnmasked.ncbi.bb" ||
   -s "bbi/${accAsm}.cpgIslandExt.ncbi.bb" ) );

my $groupName = $wrkDir;
my $asmType = "genbank";
$asmType = "refseq" if ( $groupName =~ m#/refseq/#);
$groupName =~ s#.*/$asmType/##;
$groupName =~ s#/.*##;
my $selfUrl = "${groupName}/${accAsm}";

# printf STDERR "%s\n", $accAsm;
my $htmlFile="${accAsm}.cpgIslands.html";
open (FH, ">$htmlFile") or die "can not write to $htmlFile";

my $em = "<em>";
my $noEm = "</em>";
my $assemblyDate = `grep -v "^#" $accAsm.names.tab | cut -f9`;
chomp $assemblyDate;
my $ncbiAssemblyId = `grep -v "^#" $accAsm.names.tab | cut -f10`;
chomp $ncbiAssemblyId;
my $agpFiles = `ls ${accAsm}.*.agp.ncbi.gz | xargs echo`;
chomp $agpFiles;
my $organism = `grep -v "^#" $accAsm.names.tab | cut -f5`;
chomp $organism;


print FH <<_EOF_
<h2>Description</h2>

<p>CpG islands are associated with genes, particularly housekeeping
genes, in vertebrates.  CpG islands are typically common near
transcription start sites and may be associated with promoter
regions.  Normally a C (cytosine) base followed immediately by a 
G (guanine) base (a CpG) is rare in
vertebrate DNA because the Cs in such an arrangement tend to be
methylated.  This methylation helps distinguish the newly synthesized
DNA strand from the parent strand, which aids in the final stages of
DNA proofreading after duplication.  However, over evolutionary time,
methylated Cs tend to turn into Ts because of spontaneous
deamination.  The result is that CpGs are relatively rare unless
there is selective pressure to keep them or a region is not methylated
for some other reason, perhaps having to do with the regulation of gene
expression.  CpG islands are regions where CpGs are present at
significantly higher levels than is typical for the genome as a whole.</p>

<p>
The unmasked version of the track displays potential CpG islands
that exist in repeat regions and would otherwise not be visible
in the repeat masked version.
</p>

<p>
By default, only the masked version of the track is displayed.  To view the
unmasked version, change the visibility settings in the track controls at
the top of this page.
</p>

<h2>Methods</h2>

<p>CpG islands were predicted by searching the sequence one base at a
time, scoring each dinucleotide (+17 for CG and -1 for others) and
identifying maximally scoring segments.  Each segment was then
evaluated for the following criteria:

<ul>
	<li>GC content of 50% or greater</li>
	<li>length greater than 200 bp</li>
	<li>ratio greater than 0.6 of observed number of CG dinucleotides to the expected number on the 
	basis of the number of Gs and Cs in the segment</li>
</ul>
</p>
<p>
The entire genome sequence, masking areas included, was
used for the construction of the  track <em>Unmasked CpG</em>.
The track <em>CpG Islands</em> is constructed on the sequence after
all masked sequence is removed.
</p>

<p>The CpG count is the number of CG dinucleotides in the island.  
The Percentage CpG is the ratio of CpG nucleotide bases
(twice the CpG count) to the length.  The ratio of observed to expected 
CpG is calculated according to the formula (cited in 
Gardiner-Garden <em>et al</em>. (1987)):

<pre>    Obs/Exp CpG = Number of CpG * N / (Number of C * Number of G)</pre>

where N = length of sequence.</p>

<h2>Credits</h2>

<p>This track was generated using a modification of a program developed by G. Miklem and L. Hillier 
(unpublished).</p>

<h2>References</h2>

<p>
Gardiner-Garden M, Frommer M.
<a href="http://www.sciencedirect.com/science/article/pii/0022283687906899" target="_blank">
CpG islands in vertebrate genomes</a>.
<em>J Mol Biol</em>. 1987 Jul 20;196(2):261-82.
PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/3656447" target="_blank">3656447</a>
</p>
_EOF_
   ;

close (FH);
