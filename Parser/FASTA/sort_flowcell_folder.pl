#!/usr/bin/perl

use strict;
use warnings;


### User params
my $usage = "$0 lane\n";

my $lane         = shift or die $usage;

### Get all subfolders (AMOS batches)
my @subfolders = glob("$lane/*/*");


### Debug output
foreach my $subfolder (@subfolders) {
	print "$lane/$subfolder/reads_0.fl\n";
}

foreach my $subfolder (@subfolders) {
	system("sort -n -k1 $subfolder/reads_0.fl > $subfolder/reads_0.tmp");
	system("mv $subfolder/reads_0.tmp $subfolder/reads_0.fl");
}

exit(0);
