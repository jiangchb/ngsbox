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
#  Module: Parser::ML::qval_2_maplist.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;

my $file;

my %CMD;
GetCom();


my %SEQFOLDER = ();
my %BUSTARDFOLDER = ();
my %SEQLANES = ();
my %MAPFILES = ();

### Read in input file information
print "Reading in input file..\n";
open FILE, $file or die "Cannot open file $file\n";
my $file_flag = 0;
while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (substr($line, 0, 1) ne "#") {
		if ($file_flag == 0) {
			$SEQFOLDER{$a[0]} = $a[2];
			$BUSTARDFOLDER{$a[0]} = $a[1];
			$SEQLANES{$a[0]} = $a[3];
		} else {
			$MAPFILES{$a[0]} .= $a[1].",";
		}
	} else {
		$file_flag = 1;
	}
}
close FILE;

foreach my $fc (keys %SEQFOLDER) {

	# Read in all qval from one flowcell -- will it work?
	my %ID_2_QCAL = ();
	my @lanes = split ",", $SEQLANES{$fc};
	print "now reading in qval from $fc....\n";

	for (my $i = 0; $i < @lanes; $i++) {
		my $lane = $lanes[$i];
		print "\tLane: ", $lane, "\n";
		my $bust_folder = $BUSTARDFOLDER{$fc};
		my $seq_folder = $SEQFOLDER{$fc};
		my @seq_files = glob("$seq_folder/s_$lane*seq.txt");
		foreach my $seqfile (@seq_files) {
			open FILE, $seqfile or die "Cannot open seqfile $seqfile\n";
			# create bust file name
			my @tmp = split "/", $seqfile;
			my $ending = substr($tmp[$#tmp], 0, length($tmp[$#tmp])-8);
			my $bustfile = $bust_folder."/".$ending."_qcal.txt";

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


			open BUST, $bustfile or die "Cannot open bust file $bustfile\n";
			my $count = 0;
			while (<BUST>) {
				chomp();
				$ID_2_QCAL{$ids[$count]} = $_;
				$count++;
			}
			close BUST or die "Cannot close bust file $bustfile\n";
		}
	}

	print "parsing map.lists....\n";
	# Parse map.list 
	my $mfstring = $MAPFILES{$fc};
        my @mfs = split ",", $mfstring;
        for (my $j = 0; $j<@mfs; $j++) {
                my $mf = $mfs[$j];
                print "\tMapfile:", $mf, "\n";
		open MAPFILE, $mf or die "Cannot open file $mf\n";
		open OUTFILE, "> $mf.qval" or die "Cannot open file $mf.qval\n";
		while (<MAPFILE>) {
			chomp();
			my @tmp = split " ", $_;
			if (defined($ID_2_QCAL{$tmp[3]})) {
				print OUTFILE $_, "\t", $ID_2_QCAL{$tmp[3]}, "\n";
			} else {
				print "not defined:>", $tmp[3], "<\n";
			}
		} 	
	}

}





sub GetCom{

        my @usage = ("\nUsage: $0 --file=folder.list 
                --file\tfolder.list file
                \n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s");

        die("Please specify input file\n")                              if not defined $CMD{file};

	$file = $CMD{file};

}

exit(0);
