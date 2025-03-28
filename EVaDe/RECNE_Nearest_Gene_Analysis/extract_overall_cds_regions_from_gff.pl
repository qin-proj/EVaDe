#!/usr/bin/perl

###############################################################################
# Description: This script processes GFF files to extract the boundaries of coding sequences (CDS) for each transcript. For transcripts with multiple CDS regions, it identifies the outermost boundaries to determine the complete coding region span.
# Input: GFF format file
# Output: Tab-delimited format with:
#   - Chromosome name
#   - Start position of complete CDS region
#   - End position of complete CDS region
#   - Transcript ID
###############################################################################

use strict;
use warnings;

my %cds_regions;
my %chr;
# Open input file (GFF format) from command line argument
open(my $fh, '<', $ARGV[0]) or die;

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^#/;
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split("\t", $line);

    # Process only CDS (Coding Sequence) features
    if ($feature eq 'CDS') {
        my ($transcript_id) = $attribute =~ /transcript_name "([^"]+)"/;
	    # Store chromosome information for this transcript
        $chr{$transcript_id}=$seqname;
          # Update CDS regions
        if (exists $cds_regions{$transcript_id}) {
    		my ($prev_start, $prev_end) = @{$cds_regions{$transcript_id}};
    		my $new_start = ($start < $prev_start) ? $start : $prev_start;
    		my $new_end = ($end > $prev_end) ? $end : $prev_end;
    		$cds_regions{$transcript_id} = [$new_start, $new_end];
	} else {
    		$cds_regions{$transcript_id} = [$start, $end];
}
    }
}

close($fh);
# Output results: chromosome, start position, end position, and transcript ID
foreach my $transcript_id (keys %cds_regions) {
    my ($start, $end) = @{$cds_regions{$transcript_id}};
    print "$chr{$transcript_id}\t$start\t$end\t$transcript_id\n";
}
