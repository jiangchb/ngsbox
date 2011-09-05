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
#  Module: Parser::SHORE::Variants::variant_coverage.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------------
# Get coverage at a set of polymorphic positions from consensus
# Written by Stephan Ossowski
# --------------------------------------------------------------------------


use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $snp_file  = "";
my $cons_file = "";

my %CMD;
GetCom();


### Read reference sequence
my %snp_map = ();
open SNP, $snp_file or die "Cannot open $snp_file\n";
while( <SNP> ) {
	chomp;

	# store SNP in hash
	my @e = split("\t", $_);

	$snp_map{"$e[1]#$e[2]"} = 1;
}
close SNP;
print STDERR "Finished reading SNPs\n\n";


open CONS, $cons_file or die "Cannot open $cons_file\n";
### Parse files
while( <CONS> ) {
	chomp;
	my ($chr, $pos, $call, $cov, $a, $c, $g, $t, $d, $n) = split("\t", $_);

	# SNP pos found
	if( exists $snp_map{"$chr#$pos"} ) {
		print "$chr\t$pos\t$call\t$cov\n"
	}
}

### Finish
close CONS;

exit(0);


### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--snp     STRING     SNP file
--cons    STRING     consensus file
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "snp=s", "cons=s");

	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify consensus file\n") unless defined($CMD{cons});

	$snp_file  = $CMD{snp};
	$cons_file = $CMD{cons};
}

