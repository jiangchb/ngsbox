#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski, Korbinian Schneeberger
#Date 		05/15/07
#Version	0.9
#Function	Map a set of short sequence reads to a reference genome allowing gaps
######################################################################################

package map_splicesite_hq;

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use mapping_splicesite;

sub map_splicesite
{
	my ($mismatches, $min_read_length, $repeat_mapping, $seedlength, $dir, $suffixtree) = @_;

	mapping_splicesite::flat2fasta(0, $dir);

	system("vmatch -q $dir/reads_0.fa -d -p -l $min_read_length -h $mismatches -seedlength $seedlength -s abbrev -showdesc 30 -noscore -noidentity $suffixtree | grep '.' > $dir/map_original_0.vm");

	mapping_splicesite::clean_vmatch_hamming_distance(0, $dir, 1);
	mapping_splicesite::best_evalue(0, $dir);
	mapping_splicesite::comp(0, $dir);
	mapping_splicesite::concatenate($dir, 0);
	mapping_splicesite::add_prb($dir);
	mapping_splicesite::occurrencies($dir, $repeat_mapping);
}
1;
