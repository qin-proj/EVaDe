#!/usr/bin/perl

###############################################################################
# Description: 
# This script screens for genes closest to RECNE .
###############################################################################


use strict;
use warnings;

my %gene_distances;

open(my $fh, '<', $ARGV[0]) or die ;

#Calculate and store minimum distances
while (my $line = <$fh>) {
    chomp $line;
    my @columns = split("\t", $line);

    # Extract gene coordinates and names
    my $gene1_start = $columns[1];  # Start position of the gene to be analyzed
    my $gene1_end = $columns[2];    # End position of the gene to be analyzed;
    my $gene1_name = $columns[4];   # gene name
    my $gene2_start = $columns[9];  # Start position of RECNE
    my $gene2_end = $columns[10];   # End position of RECNE

    my $distance;
    if ($gene1_start < $gene2_start) {
        $distance = $gene2_start - $gene1_end;
    } else {
        $distance = $gene1_start - $gene2_end;
    }
    $distance = abs($distance);

	my $diff1 = abs($gene1_start - $gene2_start);
	my $diff2 = abs($gene1_start - $gene2_end);
	my $diff3 = abs($gene1_end - $gene2_start);
	my $diff4 = abs($gene1_end - $gene2_end);

	my $min_diff = $diff1;
	$min_diff = $diff2 if $diff2 < $min_diff;
	$min_diff = $diff3 if $diff3 < $min_diff;
	$min_diff = $diff4 if $diff4 < $min_diff;
	$distance=$min_diff;

    if (exists $gene_distances{$gene2_start}) {
        if ($distance < $gene_distances{$gene2_start}) {
            $gene_distances{$gene2_start} = $distance;
        }
    } else {
        $gene_distances{$gene2_start} = $distance;
    }
}

close($fh);
open(my $fh2, '<', $ARGV[0]) or die;


#Output lines with minimum distances
while (my $line = <$fh2>) {
    chomp $line;
    my @columns = split("\t", $line);

    my $gene1_start = $columns[1];
    my $gene1_end = $columns[2];
    my $gene1_name = $columns[4];
    my $gene2_start = $columns[9];
    my $gene2_end = $columns[10];
	
my $diff1 = abs($gene1_start - $gene2_start);
        my $diff2 = abs($gene1_start - $gene2_end);
        my $diff3 = abs($gene1_end - $gene2_start);
        my $diff4 = abs($gene1_end - $gene2_end);
my $min_diff = $diff1;
        $min_diff = $diff2 if $diff2 < $min_diff;
        $min_diff = $diff3 if $diff3 < $min_diff;
        $min_diff = $diff4 if $diff4 < $min_diff;
        my $distance=$min_diff;

    $distance = abs($distance);
    if($distance==$gene_distances{$gene2_start}){
    	print"$line\t$distance\n";
    }
    }
    close($fh2);

