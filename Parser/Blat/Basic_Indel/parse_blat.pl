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
#  Module: Parser::Blat::Basic_Indel::parse_blat.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

my $maplist;
my $accfolder;
my $range;
my $readlength;

## Parse command lint
my %CMD;
GetCom();

my %ANCHOR = ();
get_anchor($maplist, \%ANCHOR);

my @files = glob($accfolder."/run*/*/*/len*/*blt");

print STDERR join ",", @files;

if (not -e "$accfolder/AlignmentFolder") {
	mkdir("$accfolder/AlignmentFolder");
}
open OUT, ">".$accfolder."/AlignmentFolder/map.blt";

foreach my $file (@files) {
	print STDERR $file, "\n";
	
	open FILE, $file or die "Cannot open file\n";

	while (my $line = <FILE>) {
		my @a = split " ", $line;
		my $target_chr = $a[13];
		my $target_start = $a[15];
		my $ori = $a[8];
		if (defined($ANCHOR{$a[9]})) {
			my @b = split "#", $ANCHOR{$a[9]};
			if ($target_chr == $b[0] and $target_start >= $b[1] and $target_start <= $b[2] and (($ori eq "+" and $b[3] eq "P") or ($ori eq "-" and $b[3] eq "D"))) {
				print OUT $line;
			}
		}
	}

	close FILE;
}

##########################################################################################
## parse anchor

sub get_anchor {
	my ($file, $ref) = @_;

	open FILE, $file or die "Cannot open file: ".$file."\n";

	my $count = 0;
	while (my $line = <FILE>) {
		my @a = split " ", $line;
		$count++;
		print STDERR $a[0], "\t", $a[1], "\n" if $count%100000 == 0; 
 
		if (($a[9] == 5 || $a[9] == 8) && $a[6]==1) {

			my $value = $a[0]."#";

			if ($a[4] eq "D") {
				my $begin = $a[1] - 10;
				my $end = $a[1] + $range - $readlength + 10;
				my $ori = $a[4];
				$value .= $begin."#".$end."#".$ori;
				${$ref}{$a[3]} = $value;
			}
			else {
				my $begin = $a[1] - $range + $a[7] - 10;
                                my $end = $a[1] + 10;
				my $ori = $a[4];
                                $value .= $begin."#".$end."#".$ori;
                                ${$ref}{$a[3]} = $value;
			}
		}	
	}

}


sub GetCom {

        my @usage = ("$0 --maplist file --acc accession_folder --range size --readlength length\nrange=max clone size\nreadlengh=avg sequencing length of one read (smaller than range)\n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "maplist=s", "acc=s", "range=s", "readlength=s");

        die("Please specify blat file\n") unless defined($CMD{maplist});
        die("Please specify reference file\n") unless defined($CMD{acc});
        die("Please specify range\n") unless defined($CMD{range});
        die("Please specify read length\n") unless defined($CMD{readlength});

        $maplist = $CMD{maplist};
        $accfolder = $CMD{acc};
        $range = $CMD{range};
        $readlength = $CMD{readlength};

	if ($range < $readlength) {
		die "Range needs to be larger than read length\n";
	}

}


