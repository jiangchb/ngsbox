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
#  Module: Analysis::Assembly::Scaffolding::select_contigs.pl
#  Purpose:
#  In:
#  Out:
#



### User params
my $chr              = shift;
my $trim_len         = shift;
my $min_length       = shift;
my $min_read_per_ctg = shift;
my $max_obs_exp_cov  = shift;
my $pos_shift        = shift;
my $contig_file      = shift;
my $layout_file      = shift;

### Init containers
my %ctg_seq       = ();
my %ctg_start     = ();
my %ctg_end       = ();
my %ctg_len       = ();
my %ctg_readcount = ();
my %ctg_read_bp   = ();
my %ctg_avg_cov   = ();

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


### Validate
foreach $id (sort {$ctg_start{$a}<=>$ctg_start{$b}} keys %ctg_start) {

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
	elsif(	(($ctg_avg_cov{$id} / $overall_avg_cov) >= $max_obs_exp_cov) #||
		#(($overall_avg_cov / $ctg_avg_cov{$id}) >= $max_obs_exp_cov)
	){
		delete $ctg_seq{$id};
		delete $ctg_start{$id};
		print STDERR "BAD COVERAGE: $id\t" . $ctg_len{$id} . "\t" . $ctg_avg_cov{$id} . "\n";
	}
}


### Print validated contigs
my $total_genome_coverage = 0;
my $counter = 1;
foreach $id (sort {$ctg_start{$a}<=>$ctg_start{$b}} keys %ctg_start) {
	$ctg_start{$id} += $trim_len;
	$ctg_end{$id} -= $trim_len;
	$ctg_start{$id} += $pos_shift;
	$ctg_end{$id} += $pos_shift;
	$ctg_len{$id} -= (2 * $trim_len);
	$total_genome_coverage += $ctg_len{$id};
	my $trimmed_contig = substr($ctg_seq{$id}, $trim_len, $ctg_len{$id});

	print ">contig_$chr" . "_$counter | " . $ctg_start{$id} . " | " . $ctg_end{$id} . " | " . $ctg_len{$id} . "\n" . $trimmed_contig . "\n";
	$counter++;
}

print STDERR "Average contig coverage: $overall_avg_cov\n";
print STDERR "Total genome covered: $total_genome_coverage\n";


exit(0);
