#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Std;
use warnings;
#use Bio::DB::EUtilities;
use Bio::DB::GenBank;
use LWP::Simple;

our $opt_f;
&getopts('f:');
my $infile = $opt_f or die "Provide input fasta file (aa) with -f\n";
my $nuclfile = "nucl.fa";
my $fileformat = 'Fasta';
my $framebotfile = "framebot.fa";

# The framebot.fa file for a Xander model is a peptide fasta with taxonomic information in the header
# Example framebot header:
# >134265192_Geobacillus_thermodenitrificans_NG80-2       Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Geobacillus
# This should be a comprehensive file of full-length sequences to help correctly identify the start of the gene.
# Note that the hmm seed file will likely be less comprehensive and not include positions that vary, and thus
# should not be used as the framebot set.

# The nucl.fa file is used to determine if identified sequences are chimeric, and thus should be comprehensive and
# contain only full-length sequences.

# Read a fasta file
my $in = Bio::SeqIO->new(-file => $infile);
my @accs;
my %LOOKUP;
while (my $seqo = $in->next_seq) {
    my $display_id = $seqo->display_id;
    if ($display_id =~ /(gb|emb|sp|dbj|pir|gi)\|([\w]+)/) {
	my $id = $2;
	push (@accs, $id);
	if ($display_id =~ /ref\|([\w]+)/) {
	    my $acc = $1;
	    $LOOKUP{$id} = $acc;
	}
    } elsif ($display_id =~ /tr\|(\w+)/) {
	my $record = get("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=UniProtKB&id=$1&style=raw");
	if ($record =~ /EMBL\:(\w+)/) {
	    push @accs, $1;
	} else {

	    warn "Can't find an EMBL record for $display_id\n";
	}
    } else {
	warn "Can't yet fetch nt sequence for $display_id\n";
    }
}

# Pull the nucleotide sequence and write the nucl.fa file
# NCBI pep file header:
# gi|1008832678|gb|AMR06617.1|
# ENA REST URL example:
# http://www.ebi.ac.uk/ena/data/view/AMR06617&display=fasta
# using LWP::Simple because I don't understand REST::Client
#print STDERR "Getting nucleotide sequences from EBI\n";
open my $N, ">$nuclfile";

### This code was used to hit ENA
#my $ENA_url_base = "http://www.ebi.ac.uk/ena/data/view/";
### Now that I've switched to using nrRefSeq as my reference,
### Try to use Genbank instead
### Found some code on BioStars by davidrohanedgell
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
foreach my $id (@accs) {

### This code was used to hit ENA
#    my $seq = get($ENA_url_base . $id . "&display=fasta");
#    print $N $seq unless ($seq =~ /not found/);

### This will hit NCBI
    # use the protein_id to get the links record
    my $url = $base . "elink.fcgi?dbfrom=protein&db=nuccore&id=$id";
#    my $url = $base . "elink.fcgi?dbfrom=protein&db=genome&id=$id";
#    my $url = $base . "elink.fcgi?dbfrom=gene&db=nucleotide&id=$id&linkname=gene_nuccore_refseqrna";
    my $output = get($url);
    # from the links record, extract the gi
    my $linked_gi;
    while ($output =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/sg) {
	my $linkset = $1;
	if ($linkset =~ /<Id>(\d+)<\/Id>/sg) {
	    $linked_gi = $1;
	}
    }
    # from the gi, get the CDS nucleotide fasta
    my $url3 = $base . "efetch.fcgi?db=nuccore&id=$linked_gi&rettype=fasta_cds_na";
    my $docsum = get($url3);
    $docsum =~ /(\>[^>]*${LOOKUP{$id}}[^>]*\n)\>/;
    my $seq = $1;

    print $N $seq;
    print STDERR ".";
}
close $N;
print STDERR "\n";


# Pull the full records from Genbank
#print STDERR "Getting full Genbank records for:\n" . join(" ", @accs) . "\n";
my $gb = new Bio::DB::GenBank;
my $gbo = $gb->get_Stream_by_id(\@accs);
my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                       -db      => 'taxonomy',
#                                       -rettype => 'gb',
                                       -email   => 'william.nelson@pnnl.gov',
                                       -id      => \@accs);
$factory->get_Response(-file => ".${$}_tmp_efetch.xml");

#my $gbo = Bio::SeqIO->new(-file => ".${$}_tmp_efetch");

# write the tax info to the framebot file
print STDERR "Combining sequence and taxonomy.\n";
my $frameboto = Bio::SeqIO->new(-file => ">$framebotfile",
				-format => 'Fasta');
while (my $pro = $gbo->next_seq) {
    print STDERR ",";
    my @tax = reverse($pro->species->classification);
    my $species = $tax[$#tax];
    $species =~ s/\s+/\_/g;
    $pro->display_id($pro->display_id . "_" . $species);
    $pro->description(join(";", @tax));
    $frameboto->write_seq($pro);
}
print STDERR "\n";
exit();
