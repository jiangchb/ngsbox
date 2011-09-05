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
#  Module: Support::Misc::scp_tiles.pl
#  Purpose:
#  In:
#  Out:
#



my $lane = shift;
my $filter = shift;

mkdir("tmp_4_copy");
mkdir("tmp_4_copy/L00".$lane);

for (my $e=1; $e<=40; $e++) {
	mkdir("tmp_4_copy/L00".$lane."/C$e.1");
}

for (my $c = 1; $c<=40; $c++) {
	my $sys = "cp C".$c.".1/".$filter."* tmp_4_copy/L00".$lane."/C".$c.".1/.";
	system($sys);
}

chdir("tmp_4_copy");

`scp -r * korbinian\@cgw.tuebingen.mpg.de:~/Images/.`;

chdir("..");
#system("rm -rf tmp_4_copy");


