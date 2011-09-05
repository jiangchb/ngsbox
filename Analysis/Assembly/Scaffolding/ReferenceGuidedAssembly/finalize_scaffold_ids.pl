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
#  Module: Analysis::Assembly::Scaffolding::ReferenceGuidedAssembly::finalize_scaffold_ids.pl
#  Purpose:
#  In:
#  Out:
#

my $usage = "$0 out.gp scaffoldFile chrID offset\n";

my $out = shift or die $usage;
my $file = shift or die $usage;
my $chrid = shift or die $usage;
my $offset = shift; 
die $usage if (not defined($offset));

### MUMmer output reading
my %IDEXT = ();
my %POS = ();
my $counter_id = 1+$offset;
open FILE, $out or die "Cannot open file:".$out."\n";
while (my $line = <FILE>) {
	if (substr($line, 0, 7) eq " \"Scaff") {
		my @a = split " ", $line;
		$IDEXT{substr($a[0], 1, length($a[0])-2)} = ">Scaffold_".$counter_id." | ".$chrid;
		$POS{substr($a[0], 1, length($a[0])-2)} = substr($a[1], 0, length($a[1])-1);
                $counter_id++;
	} 
	elsif (substr($line, 3, 5) eq  "Scaff") {
		my @a = split " ", $line;
                $IDEXT{substr($a[0], 2, length($a[0])-3)} = ">Scaffold_".$counter_id." | ".$chrid;
                $POS{substr($a[0], 2, length($a[0])-3)} = substr($a[1], 0, length($a[1])-1);
                $counter_id++;
	}
}
close FILE;


### Readin scaffolds
my %SEQ= ();
my $seq = "";
my $id = "";
my $na = 1;
open FILE, $file or die "Cannot open file:".$file."\n";
while (my $line = <FILE>) {
	if (substr($line, 0, 1) eq ">") {
		if ($seq ne "") {
			if (defined($IDEXT{$id})) {
				$SEQ{$POS{$id}} .= $IDEXT{$id}."\n".$seq;
			}
			else {
				$SEQ{9999999999999} .= ">Scaffold_".$counter_id." | ".$chrid." not_anchored\n".$seq;
				$counter_id++;
			}
		}
		$seq = "";
		chomp($line);
		$id = substr($line, 1, length($line)-1);
	}
	else {
		$seq .= $line
	}
}
if ($seq ne "") {
	if (defined($IDEXT{$id})) {
		$SEQ{$POS{$id}} .= $IDEXT{$id}."\n".$seq;
	}
        else {
		$SEQ{9999999999999} .= ">Scaffold_".$counter_id." | ".$chrid." not_anchored\n".$seq;
		$counter_id++;
        }
}


foreach my $scaff (sort {$a <=> $b} keys %SEQ) {
	print $SEQ{$scaff};
}

