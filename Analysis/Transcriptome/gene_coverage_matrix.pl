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
#  Module: Analysis::Transcriptome::gene_coverage_matrix.pl
#  Purpose:
#  In:
#  Out:
#


# --------------------------------------------------------------------------
# Build a genome matrix of hundreds of transcriptomes
# Written by Stephan Ossowski
# --------------------------------------------------------------------------


use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $num_samples = "";
my $trans      = "";
my $folder_list = "";

my %CMD;
GetCom();

my %gene_cov = ();
my %gene_len = ();


### Read reference sequence
open TRANS, $trans or die "Cannot open $trans\n";

# variables
my $last_gene = "";
my $last_offset = 0;

# store first entry
my $first_line = <TRANS>;
my ($id, $gene, $junk, $offset) = split(" ", $first_line);
$last_gene = $gene;
$last_offset = $offset;

# process each entry
while( <TRANS> ) {
        chomp;
	($id, $gene, $junk, $offset) = split(" ", $_);

	my @tmp = ();
	for(my $x = 0; $x < $num_samples; $x++) {
		push( @tmp, 0 );
	}

	$gene_cov{$last_gene} = \@tmp;
	$gene_len{$last_gene} = $offset - $last_offset;

	$last_gene = $gene;
	$last_offset = $offset;
}

# store last entry
my @tmp = ();
for(my $x = 0; $x < $num_samples; $x++) {
	push( @tmp, 0 );
}
$gene_cov{$last_gene} = \@tmp;
$gene_len{$last_gene} = 1000;

close TRANS;



### Prepare database tables and ecotype list
my $create_table = "CREATE table gene_coverage_matrix( gene varchar(20), gene_length int(5)";
my $eco_list = "";
my $counter = 0;


### Read all result files from project folders
open FOLDER, $folder_list or die "Cannot open $folder_list\n\n";
while( <FOLDER> ) {
	chomp;
	my $folder_name = $_;

	### Open consensus file
	open CONS, "$folder_name/consensus_summary.txt" or next;


	### Prepare database tables and ecotype list
	my @folder_structure = split("\/", $folder_name);
	my $eco = $folder_structure[0];
	$eco_list .= ($counter + 3) . "\t$eco\n";
	$create_table .= ", " . $eco . "_cov int(6)";


	### Parse files and store coverage (all and nonrep)
	while( <CONS> ) {
		chomp;
		my ($gene, $pos, $ref, $cov) = split("\t", $_);
		$gene_cov{$gene}[$counter] += $cov;
	}

	### Finish processing file
	close CONS;
	$counter++;
}


### Print coverage frequency matrix
open OUT, ">gene_coverage_matrix.txt" or die "Cannot open output file\n";
foreach my $gene (sort keys %gene_cov) {
	print OUT "$gene\t" . $gene_len{$gene};
	for(my $i = 0; $i < $counter; $i++) {
		print OUT "\t" . $gene_cov{$gene}[$i];
	}
	print OUT "\n";
}


### Print database create statement
$create_table .= ", PRIMARY KEY(gene));";
open CREATE, ">gene_coverage_matrix_create_table.txt" or die "Cannot open create file\n";
print CREATE "$create_table";
close CREATE;


### Print list of ecotypes
open LIST, ">gene_coverage_matrix_ecotype_list.txt" or die "Cannot open ectype file\n";
print LIST $eco_list;
close LIST;


### Finish
exit(0);



### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--samples    INT        Number of samples t be processed
--trans      STRING     Reference sequence trans file
--folder     STRING     Text file with folder list of resequencing projects
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "samples=s", "trans=s", "folder=s");

	die("Please specify number of samples\n") unless defined($CMD{samples});
        die("Please specify ref-seq file\n") unless defined($CMD{trans});
	die("Please specify folder file\n") unless defined($CMD{folder});

	$num_samples = $CMD{samples};
	$trans       = $CMD{trans};
	$folder_list = $CMD{folder};
}
