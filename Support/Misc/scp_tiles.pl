#! /usr/bin/perl

use strict;
use warnings;

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


