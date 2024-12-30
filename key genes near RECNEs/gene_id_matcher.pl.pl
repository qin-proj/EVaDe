#!/usr/bin/perl -w
use strict;
open IN1,"$ARGV[0]" or die;
open IN2,"$ARGV[1]" or die;
open OUT1,">$ARGV[2]" or die;
open OUT2,">$ARGV[3]" or die;

my%hash;
while(<IN1>){
	chomp;
	my@arr=split/\t/;
	my ($gene_type) = $arr[8] =~ /gene_type "([^"]+)"/;
	my ($transcript_id) = $arr[8] =~ /gene_id "([^"]+)\./;
	my ($transcript_name) = $arr[8] =~ /gene_name "([^"]+)"/;
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
