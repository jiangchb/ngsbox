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
#  Module: Parser::SAM::SAM_Indel.pl
#  Purpose:
#  In:
#  Out:
#



my $file = shift or die "\n\nUsage: $0 SAMfile\n\n";
open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	my $is_paired = 1;

	if ( $a[1] & hex(0x0004) ) {
		$is_paired = 0;
	}
	if ( $a[1] & hex(0x0008) ) {
		$is_paired = 0;
	}
	if( $a[8] == 0 ) {
		$is_paired = 0;
	}

	if($is_paired == 1) {
		print $_;
	}
}

close FILE;

exit(0);
