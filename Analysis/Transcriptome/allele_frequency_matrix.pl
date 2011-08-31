#!/usr/bin/perl

# --------------------------------------------------------------------------
# Build a genome matrix of hundreds of resequenced genomes
# Written by Stephan Ossowski
# --------------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $num_samples = "";
my $snp_file    = "";
my $folder_list = "";

my %CMD;
GetCom();

### Read reference sequence
my %snp_map = ();
open SNP, $snp_file or die "Cannot open $snp_file\n";
while( <SNP> ) {
	chomp;

	# store SNP in hash
	my @e = split("\t", $_);
	my @tmp = ();
	
	push(@tmp, $e[1]); push(@tmp, $e[2]); push(@tmp, $e[3]); push(@tmp, $e[4]); push(@tmp, $e[5]); push(@tmp, $e[6]); push(@tmp, $e[7]);

	for(my $i = 0; $i < $num_samples; $i++) {
		push(@tmp, 0);
		push(@tmp, 0);
	}
	
	$snp_map{"$e[1]#$e[2]"} = \@tmp;
}
close SNP;
print STDERR "Finished reading SNPs\n\n";


### Prepare database tables and ecotype list
my $create_table = "CREATE table allele_frequency_matrix( chromosome tinyint(1), position int(9), ref char(1), variant char(1), quality tinyint(2), support int(6), concordance double";
my $eco_list = "";
my $counter = 7;


### Read all result files from project folders
open FOLDER, $folder_list or die "Cannot open $folder_list\n\n";
while( <FOLDER> ) {
	chomp;
	my $folder_name = $_;

	### Open files
	open CONS, "$folder_name/consensus_summary.txt" or next;

	### Prepare database tables and ecotype list
	my @folder_structure = split("\/", $folder_name);
	my $eco = $folder_structure[0];
	$eco_list .= ($counter + 1) . "\t$eco reference count\n";
	$eco_list .= ($counter + 2) . "\t$eco variant count\n";
	$create_table .= ", " . $eco . "_ref int(6), " . $eco . "_var int(6)";
	print STDERR "Now starting to process sample $eco\n";

	### Initialize variables
	my $chr = "";
	my $pos = -1;
	my $call = "";
	my $cov = -1;

	### Parse files
	while( <CONS> ) {
		my $line = $_;
		chomp($line);

		my %allele = ();	
		($chr, $pos, $call, $cov, $allele{'A'}, $allele{'C'}, $allele{'G'}, $allele{'T'}, $allele{'-'}) = split("\t", $line);

		# Reference call found
		if( exists $snp_map{"$chr#$pos"} ) {
			my $ref = $snp_map{"$chr#$pos"}[2];
			my $var = $snp_map{"$chr#$pos"}[3];

			$snp_map{"$chr#$pos"}[$counter]   = $allele{$ref};
			$snp_map{"$chr#$pos"}[$counter+1] = $allele{$var};
		}
	}

	### Finish
	close CONS;
	$counter += 2;
}

### Print allele frequency matrix
open OUT, ">allele_frequency_matrix.txt" or die "Cannot open output file\n";
foreach my $chrpos (sort keys %snp_map) {
	for(my $i = 0; $i < $counter; $i++) {
		print OUT $snp_map{$chrpos}[$i] . "\t";
	}
	print OUT "\n";
}

### Print database create statement
$create_table .= ", PRIMARY KEY(chromosome, position));";
open CREATE, ">allele_frequency_matrix_create_table.txt" or die "Cannot open create file\n";
print CREATE "$create_table";
close CREATE;

### Print list of ecotypes
open LIST, ">allele_frequency_matrix_ecotype_list.txt" or die "Cannot open ectype file\n";
print LIST $eco_list;
close LIST;

exit(0);


### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--samples    INT        Number of samples t be processed
--snp        STRING     SNP file
--folder     STRING     Text file with folder list of resequencing projects
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "samples=s", "snp=s", "folder=s");

	die("Please specify number of samples\n") unless defined($CMD{samples});
	die("Please specify snp file\n") unless defined($CMD{snp});
	die("Please specify folder file\n") unless defined($CMD{folder});

	$num_samples = $CMD{samples};
	$snp_file    = $CMD{snp};
	$folder_list = $CMD{folder};
}

