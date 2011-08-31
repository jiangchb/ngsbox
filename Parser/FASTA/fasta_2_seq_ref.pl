#!/usr/bin/perl -W
# written by korbinian schneeberger


use warnings;
use strict;
use Getopt::Long;

my %CMD;
my $assem;
my $file;
my $chr;
my $chr_count = 0;
my $pos = 1;

GetCom();

open FILE, $CMD{file} or die "Cannot open file:", $CMD{file},"\n";
open REF, ">".$CMD{file}.".seq_ref";
open MAX, ">".$CMD{file}.".seq_max";
open FAS, ">".$CMD{file}.".fa";

RUN: while (my $l = <FILE>) {
	chomp($l);
	if (substr($l, 0, 1) eq ">") {
		if (defined($chr)) {
			print MAX $assem, "\t", $chr_count, "\t", $chr, "\t", "\\N", "\t", $pos, "\n"; 
		}
		my @a = split " ", $l;
		$chr = substr($a[0], 1, length($a[0])-1);
		$chr_count++;
		$pos = 1;
		print FAS ">$chr_count\n";
	} else {
		if (defined($chr)) {
			for (my $i=0; $i<length($l); $i++) {
				print REF $chr_count, "\t", $pos, "\t", substr($l, $i, 1), "\n";
				$pos++;
			}
			print FAS "$l\n";
		}
	}
}
if (defined($chr)) {
	print MAX $assem, "\t", $chr_count, "\t", $chr, "\t", "\\N", "\t", $pos, "\n";
}

close FILE; close REF; close MAX; close FAS;
exit(0);

sub GetCom{

	my @usage = ("\nUsage: $0
                --file=file\treference genome
		--assem=name\tassembly name

	  	description:
		Makes a seq_ref and seq_max table from a fasta file
		\n\n");
 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file=s", "assem=s");

	die("Please specify file\n")	if not defined $CMD{file};
	die("Please specify assem\n")    if not defined $CMD{assem};

	$file = $CMD{file};
	$assem = $CMD{assem};

}
	      

