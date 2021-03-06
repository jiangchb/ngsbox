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
#  Module: Analysis::Eighties::genome_matrix::creation::genome_matrix.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------------
# Build a genome matrix of hundreds of resequenced genomes
# Written by Stephan Ossowski
# --------------------------------------------------------------------------


use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $refseq      = "";
my $folder_list = "";

my %CMD;
GetCom();

### Open output file
open OUT, ">genome_matrix.txt" or die "Cannot open output file\n";

### Read reference sequence
my $chr = 1;
my $seq = "";
open REFSEQ, $refseq or die "Cannot open $refseq\n";
while( <REFSEQ> ) {
	chomp;
	if($_ =~ /^>/ ) {
		if($seq ne "") {
			for(my $i = 1; $i <= length($seq); $i++) {
				print OUT "$chr\t$i\t" . substr($seq, $i-1, 1) . "\n";
			}
		}

		$chr = $_;
		$chr =~ s/>//g;
		$seq = "";
	}
	else {
		$seq .= $_;
	}
}
if($seq ne "") {
	for(my $i = 1; $i <= length($seq); $i++) {
		print OUT "$chr\t$i\t" . substr($seq, $i-1, 1) . "\n";
	}
}

close REFSEQ;
close OUT;

### Prepare database tables and ecotype list
my $create_table = "CREATE table genome_matrix( chromosome tinyint(1), position int(9), ref char(1)";
my $eco_list = "";
my $counter = 1;

### Read all result files from project folders
open FOLDER, $folder_list or die "Cannot open $folder_list\n\n";
while( <FOLDER> ) {
	chomp;
	my $folder_name = $_;

	### Open files
	open RCALL, "$folder_name/quality_reference.txt" or next;
	open VAR, "$folder_name/quality_variant.txt" or die "Cannot open variant calls\n";
	open DONE, "genome_matrix.txt" or die "Cannot open done file\n";
	open OUT, ">genome_matrix2.txt" or die "Cannot open outfile\n";

	### Prepare database tables and ecotype list
	my @folder_structure = split("\/", $folder_name);
	my $eco = $folder_structure[0];
	$eco_list .= "$counter\t$eco\t";
	$eco =~ s/-/_/g;
	$eco_list .= "$eco\n";
	$counter++;
	$create_table .= ", " . $eco . "_call char(1), " . $eco . "_qual tinyint(2)";

	### Initialize variables
	$chr = 1;
	my $current_pos = 1;
	my $rcall_line = <RCALL>;
	chomp($rcall_line);
	my $var_line = <VAR>;
	chomp($var_line);

	### Parse files
	while( <DONE> ) {
		my $done_line = $_;
		chomp($done_line);

		my ($done1, $done2) = split("\t", $done_line);
		my @rcall = split("\t", $rcall_line);
		my @var = split("\t", $var_line);

		my $quality = 0;
		my $call = "N";

		# Reference call found
		if( ($done1 eq $rcall[1]) && ($done2 eq $rcall[2]) ) {
			$quality = $rcall[5];
			$call = $rcall[3];
			$rcall_line = <RCALL>;
			if(! defined $rcall_line) { $rcall_line = "-1\t-1\t-1"; }
			chomp($rcall_line);
		}

		# Variant found (SNP or deletion)
		if( ($done1 eq $var[1]) && ($done2 eq $var[2]) ) {
			$quality = $var[5];
			$call = $var[4];
			$var_line = <VAR>;
			if(! defined $var_line) { $var_line = "-1\t-1\t-1"; }
			chomp($var_line);
		}

		### Print position
		print OUT "$done_line\t$call\t$quality\n";
	}

	### Close all files
	close RCALL;
	close VAR;
	close DONE;
	close OUT;

	### mv temporary output file to main output file
	system("mv genome_matrix2.txt genome_matrix.txt");
}

### Print database create statement
$create_table .= ", PRIMARY KEY(chromosome, position));";
open CREATE, ">genome_matrix_create_table.txt" or die "Cannot open create file\n";
print CREATE "$create_table";
close CREATE;

### Print list of ecotypes
open LIST, ">genome_matrix_ecotype_list.txt" or die "Cannot open ectype file\n";
print LIST $eco_list;
close LIST;

exit(0);

### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--refseq     STRING     Reference sequence file
--folder     STRING     Text file with folder list of resequencing projects
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "refseq=s", "folder=s");

	die("Please specify refseq file\n") unless defined($CMD{refseq});
	die("Please specify folder file\n") unless defined($CMD{folder});

	$refseq      = $CMD{refseq};
	$folder_list = $CMD{folder};
}

