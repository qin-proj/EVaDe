#!/usr/bin/perl -w
###############################################################################
# Description: Matches gene IDs between different annotation versions.
# The script matches records based on chromosome and gene ID, and separates matched and unmatched entries into two output files.
###############################################################################


use strict;
open IN1,"$ARGV[0]" or die; # Input GFF file
open IN2,"$ARGV[1]" or die;	# Input coordinates file
open OUT1,">$ARGV[2]" or die;	# Output file for matched entries
open OUT2,">$ARGV[3]" or die;	# Output file for unmatched entries


my%hash;
while(<IN1>){
	chomp;
	my@arr=split/\t/;

	# Extract gene information from the attribute field (column 9)
	my ($gene_type) = $arr[8] =~ /gene_type "([^"]+)"/;
	my ($transcript_id) = $arr[8] =~ /gene_id "([^"]+)\./;
	my ($transcript_name) = $arr[8] =~ /gene_name "([^"]+)"/;
	 # Create position string with chromosome, start, end, gene type, and name
	my $pos="$arr[0]\t$arr[3]\t$arr[4]\t$gene_type\t$transcript_name";
	my$name="$arr[0]-$transcript_id";
	$hash{$name}=$pos;
}
while(<IN2>){
	chomp;
	my@arr=split/\t/;
	my$nam="$arr[0]-$arr[5]";
	if (defined $hash{$nam}){
			print OUT1"$_\t$hash{$nam}\n";
		}else{
			print OUT2"$_\n";
			}

}
close IN1;
close IN2;
close OUT1;
close OUT2;
