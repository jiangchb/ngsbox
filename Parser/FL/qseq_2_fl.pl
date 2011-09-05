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
#  Module: Parser::FL::qseq_2_fl.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;

my $bustardfolder;
my $geraldfolder;
my $file;

my %CMD;
GetCom();

my @IDS = ();
my %ID_2_QCAL = ();



my @seq_files = glob("$bustardfolder/*seq.txt");
foreach my $seqfile (@seq_files) {
	open FILE, $seqfile or die "Cannot open seqfile $seqfile\n";
	# create gerald file name
	my @tmp = split "/", $seqfile;
	my $ending = substr($tmp[$#tmp], 0, length($tmp[$#tmp])-8);
	my $geraldfile = $geraldfolder."/".$ending."_qcal.txt";

print STDERR $ending, "\n";

	my @ids = ();
	while (my $line = <FILE>) {
		my @a = split " ", $line;
		my $id = $a[0];
		for (my $i=1; $i<=3; $i++) { 
			if ($a[$i] eq "-0") {
				$a[$i] = "0";
			}
			$id .= "000" if length($a[$i]) == 1;
			$id .= "00" if length($a[$i]) == 2;
			$id .= "0" if length($a[$i]) == 3;
			$id .= $a[$i];
		}
		push @ids, $id;
	}
	close FILE or die "Cannot close seqfile \n";


	open GERALD, $geraldfile or die "Cannot open geraldfile $geraldfile\n";
	my $count = 0;
	while (<GERALD>) {
		chomp();
		$ID_2_QCAL{$ids[$count]} = $_;
		$count++;
	}
	close GERALD or die "Cannot close bust file $geraldfile\n";
}

print "parsing fl file....\n";

open FLFILE, $file or die "Cannot open file $file\n";
while (<FLFILE>) {
	chomp();
	my @tmp = split " ", $_;
	if (defined($ID_2_QCAL{$tmp[0]})) {
		print $_, "\t", $ID_2_QCAL{$tmp[0]}, "\n";
	} else {
		print STDERR "not defined:>", $tmp[0], "<\n";
	}
} 	





sub GetCom{

        my @usage = ("\nUsage: $0 --file=fl-file --geraldfolder=directoy --bustardfolder=directory
		--geraldfolder\t
		--bustardfolder\t
                --file\tfolder.list file

		!! Prints to stdout !!

                \n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "geraldfolder=s", "bustardfolder=s");

        die("Please specify input file\n")                              if not defined $CMD{file};
	die("Please specify geraldfolder\n")                              if not defined $CMD{geraldfolder};
	die("Please specify bustardfolder\n")                              if not defined $CMD{bustardfolder};

	$file = $CMD{file};
	$geraldfolder = $CMD{geraldfolder};
	$bustardfolder = $CMD{bustardfolder};

}

exit(0);
