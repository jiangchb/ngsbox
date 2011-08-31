#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Add prb, qval and chas to map.list
######################################################################################

use strict;
use warnings;

my $dir = shift;

system("cut -f 1-9 $dir/map.list > $dir/map.cut");
system("sort --buffer-size=30% -n -k4 $dir/map.cut > $dir/map.resort");

open (READS, "$dir/reads_0.fl") or die "cannot open READS\n";
open (MAP, "$dir/map.resort") or die "cannot open MAP\n";
open (OUT, ">$dir/map.list.new.unsorted") or die "cannot open OUT\n";


my $map_line = <MAP>;
chomp($map_line);
my @map_elem = split(/\t/, $map_line);


### Loop through original reads and get prb entry
while( my $read_line = <READS> ) {
	chomp($read_line);
	my ($read_id, $read_seq, $prb, $qCal, $chas) = split("\t", $read_line);

	if( $map_elem[3] eq $read_id) {
		while($map_elem[3] eq $read_id) {
			print OUT "$map_line\t$prb\t$qCal\t$chas\n";
			$map_line = <MAP>;
			if(not defined $map_line) { $map_elem[3] = "endoffile"; last; }
			chomp($map_line);
			@map_elem = split(/\t/, $map_line);
		}
	}
}

close READS; close MAP; close OUT;

system("sort --buffer-size=30% -n -k1 -k2 $dir/map.list.new.unsorted > $dir/map.list.new");
system("rm $dir/map.cut $dir/map.resort $dir/map.list.new.unsorted");

exit(0)
