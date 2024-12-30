#!/usr/bin/perl

use strict;
use warnings;

my %cds_regions;
my %chr;
open(my $fh, '<', $ARGV[0]) or die;

while (my $line = <$fh>) {
    chomp $line;
    next if $line =~ /^#/;
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split("\t", $line);
    if ($feature eq 'CDS') {
        my ($transcript_id) = $attribute =~ /transcript_name "([^"]+)"/;
	$chr{$transcript_id}=$seqname;
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

foreach my $transcript_id (keys %cds_regions) {
    my ($start, $end) = @{$cds_regions{$transcript_id}};
    print "$chr{$transcript_id}\t$start\t$end\t$transcript_id\n";
}
