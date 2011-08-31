#!/usr/bin/perl 

our $VERSION = '1.0';

use strict;
use warnings;

use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

use random_snp_set;
use random_deletion_set;
use random_insertion_set;


### User parameters
my $reference_file = "";	# Reference chromosome
my $deletion_file  = "";	# Simulated deletions
my $insertion_file = "";	# Simulated insertions
my $snp_file       = "";	# Simulated SNPs
my $out_file       = "";	# Outfile for simulated chromosome

my %CMD;
GetCom();


### Read sets of SNPs and indels
my $snp_creator = new random_snp_set();
my %snp = $snp_creator->read("$snp_file");
my $del_creator = new random_deletion_set();
my %del = $del_creator->read("$deletion_file");
my $ins_creator = new  random_insertion_set();
my %ins = $ins_creator->read("$insertion_file");


### Read deletions by single position
my %deletion = ();
foreach my $beg( sort {$a<=>$b} keys %del ) {
	for(my $i = $beg; $i <= $del{$beg}; $i++) {
		$deletion{$i} = 1;
	}
}


### Read reference chromosome sequence
my $ref_seq = "";
open REF, $reference_file or die "Cannot open input file $reference_file\n";
while( <REF> ) {
	chomp($_);
	if($_ !~ />/) {	
		$ref_seq .= $_;
	}
}

open OUT, ">$out_file" or die "Cannot open output file $out_file\n";
print OUT ">Chr1 Left Arm\n";

# For each position in the original chromosome
for(my $i = 0; $i < length($ref_seq); $i++) {

	my $pos = $i + 1;
	my $base = substr($ref_seq, $i, 1);

	if($pos == 13700000) { print OUT "\n>Chr1 Left Arm\n"; }
	if( ($pos >= 13700000) && ($pos <= 15900000) ) { next; }
	
	# Case insertion
	if(exists $ins{$pos}) {
		print OUT $ins{$pos};
	}

	# IF position is not deleted
	if(! exists $deletion{$pos}) {

		# Case SNP
		if(exists $snp{$pos}) {
			$base = $snp{$pos};
		}

		# Print current base
		print OUT "$base";

		if($pos % 60 == 0) { print OUT "\n"; }
	}
}

exit(0);


### Read command line parameters
sub GetCom {

	my @usage = ("$0

Mandatory:
--reference    STRING    Reference chromosome
--delfile      STRING    Infile for deletions
--insfile      STRING    Infile for insertions
--snpfile      STRING    Infile for SNPs
--outfile      STRING    Outfile for simulated chromosome

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "reference=s", "delfile=s", "insfile=s", "snpfile=s", "outfile=s");
	
	die("Please specify reference infile\n") unless defined($CMD{reference});
	die("Please specify deletion infile\n") unless defined($CMD{delfile});
	die("Please specify insertion infile\n") unless defined($CMD{insfile});
	die("Please specify SNP infile\n") unless defined($CMD{snpfile});
	die("Please specify outfile\n") unless defined($CMD{outfile});

	$reference_file  = $CMD{reference};
	$deletion_file   = $CMD{delfile};
	$insertion_file  = $CMD{insfile};
	$snp_file        = $CMD{snpfile};
	$out_file        = $CMD{outfile};
}

