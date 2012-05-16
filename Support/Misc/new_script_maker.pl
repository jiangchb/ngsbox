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
#  Module: Support::Misc::new_script_maker.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "$0 raw\|fasta\n";
my $type = shift;
if ($type ne "raw" and $type ne "fasta")  { die $usage; }

### print HEADER

if (1) {

print STDOUT 

"my \$usage = \"\$0 file\\n\";

open FILE, shift or die \$usage;
";

}

## print RAW

if ($type eq "raw") {

print STDOUT
"while (<FILE>) {
	my \@a = split \" \";
	
}
";

}

## print FASTA

if ($type eq "fasta") {

print STDOUT

"my \%SEQ = ();
my \$seq = \"\";
my \$id = \"\";

while (<FILE>) {
	if (substr(\$_, 0, 1) eq \">\") {
		my \@a = split \" \", \$_;
		if (\$seq ne \"\") {
			\$SEQ{\$id} = \$seq;
		}
		\$id = substr(\$a[0], 1, length(\$a[0])-1);
		\$seq = \"\";
	}
	else {
		chomp(\$_);
		\$seq .= \$_;
	}
}
if (\$seq ne \"\") {
	\$SEQ{\$id} = \$seq;
}
";

}



