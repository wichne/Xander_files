#!/usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Std;
use warnings;
use Bio::Seq;
#use LWP::Simple;

our $opt_f;
&getopts('f:');
my $infile = $opt_f or die "Provide input fasta file (aa) with -f\n";

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
my $nuclfile = Bio::SeqIO->new(-file => ">nucl.fa",
			       -format => 'Fasta');
my $framebotfile = Bio::SeqIO->new(-file => ">framebot.fa",
				   -format => 'Fasta');

while (my $seqo = $in->next_seq) {
    my $display_id = $seqo->display_id;
    my $id;
    my $acc;
    if ($display_id =~ /(gb|emb|sp|dbj|pir|gi)\|([\w]+)/) {
	$id = $2;
	if ($display_id =~ /ref\|([\w]+)/) {
	    $acc = $1;
	}
    } elsif ($display_id =~ /tr\|(\w+)/) {
	my $record = get("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=UniProtKB&id=$1&style=raw");
	if ($record =~ /EMBL\:(\w+)/) {
	    $acc = $id = $1;
	} else {
	    warn "Can't find an EMBL record for $display_id\n";
	}
    } else {
	warn "Can't yet fetch nt sequence for $display_id\n";
    }

    # for the nucl.fa file
    my $sequence = &fetch_sequence($id, $acc);
    if (! $sequence) { warn "No seq for $display_id\n"; }
    else {
	my $nuclo = Bio::Seq->new(-id    => $id,
				  -seq   => $sequence);
	$nuclfile->write_seq($nuclo);
    }

    my $lineage  = &fetch_lineage($id);
    if (! $lineage) { warn "No linage for $display_id\n"; }
    else {
	$seqo->desc($lineage);
	$framebotfile->write_seq($seqo);
    }
}

sub fetch_lineage {
    my $id = shift;
    # Getting taxonomy
    my $NCBI_url_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    # use the protein_id to link to nucleotide
    my $url = $NCBI_url_base . "elink.fcgi?dbfrom=protein&db=taxonomy&id=$id";
    my $output = `wget -q -O - "$url"`;
    # from the links record, extract the taxonid
    if ($output) {
	my $linked_taxonid;
	while ($output =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/sg) {
	    my $linkset = $1;
	    if ($linkset =~ /<Id>(\d+)<\/Id>/sg) {
		$linked_taxonid = $1;
		# from the taxonid, get the lineage
		my $url3 = $NCBI_url_base . "efetch.fcgi?db=taxonomy&id=$linked_taxonid";
		my $docsum = `wget -q -O - "$url3"`;
		if ($docsum =~ /<Lineage>(.*?)\n*<\/Lineage>/sg) {
		    my $lineage = $1;
		    $lineage =~ s/\;\s+/\;/g;
		    $lineage =~ s/\s+/\_/g;
		    return($lineage);
		}
	    }
	}
    }
}

sub fetch_sequence {
    my $id = shift;
    my $acc = shift;
    my $NCBI_url_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $NCBI_url_base . "elink.fcgi?dbfrom=protein&db=nuccore&id=$id";
    my $output = `wget -q -O - "$url"`;
    # from the links record, extract the gi
    if ($output) {
	my $linked_gi;
	while ($output =~ /<LinkSetDb>(.*?)<\/LinkSetDb>/sg) {
	    my $linkset = $1;
	    if ($linkset =~ /<Id>(\d+)<\/Id>/sg) {
		$linked_gi = $1;
		# from the gi, get the CDS nucleotide fasta
		my $url3 = $NCBI_url_base . "efetch.fcgi?db=nuccore&id=$linked_gi&rettype=fasta_cds_na";
		my $docsum = `wget -q -O - "$url3"`;
		if ($docsum =~ /(\>[^>]*$acc[^>]*\n)\>/) {
		    my $seq = $1;
		    $seq =~ s/^[^\n]*\n//;
		    $seq =~ s/\n//g;
		    return($seq);
		}
	    }
	}
   }
}
exit();
