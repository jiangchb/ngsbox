#!/usr/bin/perl

use strict;
use warnings;

### User params
my $trim_len         = shift;
my $min_length       = shift;
my $min_read_per_ctg = shift;
my $min_bp_per_mm    = shift;
my $max_obs_exp_cov  = shift;

my $contig_file    = shift;
my $layout_file    = shift;
my $unseq_file     = shift;
my $snp_file       = shift;
my $deletion_file  = shift;
my $insertion_file = shift;


### Init containers
my %ctg_seq       = ();
my %ctg_start     = ();
my %ctg_end       = ();
my %ctg_len       = ();
my %ctg_readcount = ();
my %ctg_read_bp   = ();
my %ctg_avg_cov   = ();

my %splitter = ();
my %unseq_N  = ();
my %snp      = ();

my $id = -1;


### Parse contigs
open CONTIG, $contig_file or die "Cannot open $contig_file\n";

while(<CONTIG>) {
	chomp($_);

	if (substr($_, 0, 1) eq ">") {
		$id = substr($_, 1);
	}
	else {
		$ctg_seq{$id} .= $_;
	}
}
close CONTIG;


### Parse contig layout file
$id = -1;
open LAYOUT, $layout_file or die "Cannot open $layout_file\n";

while( <LAYOUT> ) {
	chomp($_);
	my $line = $_;

	if($line ne "") {
		### New contig
		if( substr($line, 0, 1) eq "C")  {
			my @a = split(" ", $line);
			$id = $a[1];
			$ctg_readcount{$id} = $a[2];
			$ctg_read_bp{$id} = 0;

			# Get start and end of the contig relative to the reference genome
			if( substr($a[4], 0, 1) ne "-")  {
				my ($start, $end) = split("-", $a[4]);
				$ctg_start{$id} = $start;
				$ctg_end{$id} = $end;
			}
			else {
				my ($junk, $start, $end) = split("-", $a[4]);
				$ctg_start{$id} = -$start;
				$ctg_end{$id} = $end;
			}
		}

		### Reads contained in current contig
		else {
			my @a = split(" ", $line);
			$ctg_read_bp{$id} += (abs($a[1] - $a[2]) + 1);
		}
	}
}
close LAYOUT;


### Set contig length and avg coverage
my $total_bp  = 0;
my $total_len = 0;
foreach $id (keys %ctg_seq) {

	$ctg_len{$id} = length($ctg_seq{$id});
	$ctg_avg_cov{$id} = $ctg_read_bp{$id} / $ctg_len{$id};

	$total_len += $ctg_len{$id};
	$total_bp += $ctg_read_bp{$id};
}
my $overall_avg_cov = $total_bp / $total_len;


### Parse unsequenced regions
open UNSEQ, $unseq_file or die "Cannot open $unseq_file\n";
while( <UNSEQ> ) {
	my @a = split("\t", $_);
	my $id = $a[1];

	if(! exists $unseq_N{$id} ) {
		my @tmp1 = ();
		my @tmp2 = ();
		$unseq_N{$id} = \@tmp1;
		$splitter{$id} = \@tmp2;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		my $unseq_seq = substr($ctg_seq{$id}, $a[2] - 1, $a[3] - $a[2] + 1);

		if($unseq_seq =~ /N/) {
			push( @{$unseq_N{$id}}, "$a[2]-$a[3]" );
		}
		else {
			push( @{$splitter{$id}}, "$a[2]-$a[3]" );
		}
	}
}


### Parse SNPs
open SNP, $snp_file or die "Cannot open $snp_file\n";
while( <SNP> ) {
	my @a = split("\t", $_);
	my $id = $a[1];

	if(! exists $snp{$id} ) {
		my %tmp = ();
		$snp{$id} = \%tmp;
	}

	# Do not count in trimming range
	if( ($a[2] > $trim_len) && ($a[2] < ($ctg_len{$id} - $trim_len)) ) {
		# Select high quality SNPs
		if($a[5] >= 30) {
			$snp{$id}{$a[2]} = 1;
		}
	}
}


### Parse small deletions
open DELETION, $deletion_file or die "Cannot open $deletion_file\n";
while( <DELETION> ) {
        my @a = split("\t", $_);
	my $id = $a[1];

	if(! exists $splitter{$id} ) {
		my @tmp = ();
		$splitter{$id} = \@tmp;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		push( @{$splitter{$id}}, "$a[2]-$a[3]");
	}
}


### Parse small insertions
open INSERTION, $insertion_file or die "Cannot open $insertion_file\n";
while( <INSERTION> ) {
        my @a = split("\t", $_);
	my $id = $a[1];

	if(! exists $splitter{$id} ) {
		my @tmp = ();
		$splitter{$id} = \@tmp;
	}

	if( ($a[2] > $trim_len) && ($a[3] < ($ctg_len{$id} - $trim_len)) ) {
		push( @{$splitter{$id}}, "$a[2]-$a[3]");
	}
}


# Split gapped contigs
foreach $id ( keys %splitter ) {
	my %split_locus = ();

	foreach my $locus ( @{$splitter{$id}} ) {
		my ($s, $e) = split("-", $locus);
		$s -= 20; $e += 20;
		for(my $i = $s; $i <= $e; $i++) {
			substr($ctg_seq{$id}, $i-1, 1, "N");
			push( @{$unseq_N{$id}}, "$s-$e" );
			delete $snp{$id}{$i};
		}
	}
}


### Validate
foreach $id (sort {$ctg_start{$a}<=>$ctg_start{$b}} keys %ctg_start) {

	my $remapping_mm =  scalar( keys %{$snp{$id}} );

	# Remove short contigs
	if($ctg_len{$id} < $min_length) {
		delete $ctg_seq{$id};
		delete $ctg_start{$id};
		print STDERR "BAD LENGTH: $id\t" . $ctg_len{$id} . "\n";
	}

	# Remove low read count contigs
	elsif( $ctg_readcount{$id} < $min_read_per_ctg) {
		delete $ctg_seq{$id};
		delete $ctg_start{$id};
		print STDERR "BAD COUNT: $id\t" . $ctg_len{$id} . "\t" . $ctg_readcount{$id} . "\n";
	}

	# Remove bad coverage contigs
	elsif(	(($ctg_avg_cov{$id} / $overall_avg_cov) >= $max_obs_exp_cov) ||
		(($overall_avg_cov / $ctg_avg_cov{$id}) >= $max_obs_exp_cov)
	){
		delete $ctg_seq{$id};
		delete $ctg_start{$id};
		print STDERR "BAD COVERAGE: $id\t" . $ctg_len{$id} . "\t" . $ctg_avg_cov{$id} . "\n";
	}

	# Remove erroneous contigs
	elsif( ($remapping_mm != 0) && (($ctg_len{$id} / $remapping_mm) < $min_bp_per_mm) ) {
		delete $ctg_seq{$id};
		delete $ctg_start{$id};
		print STDERR "MISMATCHES: $id\t" . $ctg_len{$id} . "\t" . $ctg_len{$id} / $remapping_mm . "\n";
	}
}


### Print validated contigs
my $total_genome_coverage = 0;
foreach $id (sort {$ctg_start{$a}<=>$ctg_start{$b}} keys %ctg_start) {
	# TODO split at 'N' stretches

	$ctg_start{$id} += $trim_len;
	$ctg_end{$id} -= $trim_len;
	$ctg_len{$id} -= (2 * $trim_len);
	$total_genome_coverage += $ctg_len{$id};
	my $trimmed_contig = substr($ctg_seq{$id}, $trim_len, $ctg_len{$id});

	print ">$id | " . $ctg_start{$id} . " | " . $ctg_end{$id} . " | " . $ctg_len{$id} . "\n" . $trimmed_contig . "\n";
}

print STDERR "Average contig coverage: $overall_avg_cov\n";
print STDERR "Total genome covered: $total_genome_coverage\n";


exit(0);
