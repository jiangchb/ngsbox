#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Parser::ML::quality_2_maplist.pl
#  Purpose:
#  In:
#  Out:
#

######################################################################################
#Author 	Stephan Ossowski
#Date 		10/23/07
#Version	0.9
#Function	Add prb, qval and chas to map.list
######################################################################################


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
