#!/usr/bin/perl -w

use strict;

my $start = shift;
my $win   = shift;
my $file = shift;

open FILE, $file or die;
my $last_len = $start;


### Plot count in sliding window
#my $win_count = 0;
#while( <FILE> ) {
#	my ($len, $count, $sum) = split("\t", $_);
#
#	if($len >= $start) {
#		if( ($len%$win==0) && ($len > $start) ) {
#			print "$last_len\t$win_count\n";
#			$win_count = 0;
#			$last_len = $len;
#		}
#
#		$win_count+=$count;
#	}
#}


### Plot sum up to sliding window
my $new_sum = 0;
while( <FILE> ) {
	my ($len, $count, $sum) = split("\t", $_);

	if($len >= $start) {
		if( $len%$win == 0) {
			print "$len\t$new_sum\n";
		}
		$new_sum+=$count;
	}
}
