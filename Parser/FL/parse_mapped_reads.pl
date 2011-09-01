#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Read prb values
######################################################################################

use strict;
use warnings;

my $original_file = shift;
my $not_mapped_file = shift;

open (ALL, $original_file) or die "cannot open $original_file\n";
open (REST, $not_mapped_file) or die "cannot open $not_mapped_file\n";

my $rest_line = <REST>;
my ($rest_id, $rest_seq) = split("\t", $rest_line);

# Loop through original reads and get prb entry
while( my $all_line = <ALL> ) {
	chomp($all_line);
	my ($all_id, $all_seq, $all_prb, $all_qCal, $all_chas) = split(/\t/, $all_line);

	if($rest_id eq $all_id) {
		$rest_line = <REST>;
		if(! defined $rest_line) { $rest_line = "END\tNA"; }
		($rest_id, $rest_seq) = split("\t", $rest_line);
	}
	else {
		print "$all_id\t$all_seq\t$all_prb\t$all_qCal\t$all_chas\n";
	}
}

close ALL; close REST;

exit(0)
