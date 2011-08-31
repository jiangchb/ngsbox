#! /usr/bin/perl
use strict;

my $mfa = shift;

open MFA, $mfa or die "Cannot open file\n";
open CONTIG, "> original.contigs.txt";
open PRIMER3, "> primer3.contigs.txt";


#####################################
#### Parse contig ids to Primer3 format

while (my $line = <MFA>) {
	chomp();
	if ($line ne "") {
		if ($line =~ /reference/) {
			my @entries = split " ", $line;
			my $id = substr($entries[0], 5, length($entries[0])-5);
			my $col = <MFA>;
			chomp($col);
			my $junk = <MFA>;
			my $eco = <MFA>;

			# Print contig file for future comparison with Sanger sequenced stuff
			my $contig_seq = $eco;
			$contig_seq =~ s/-//g;
			print CONTIG ">$id\n$contig_seq\n";
			#print length($contig_seq), "\n";

			# Print primer3 masked sequences
			print PRIMER3 ">$id\n";
			for (my $i = 0; $i < length($eco); $i++) {
				if ($i<100 or $i > length($eco)-100) {
					if (substr($eco, $i, 1) eq substr($col, $i, 1)) {
						print PRIMER3 substr($eco, $i, 1);
					}
					else {
						print PRIMER3 "N";
					}
				}
				else {
					if (substr($eco, $i, 1) ne "-") {
						print PRIMER3 "N";
					}
				}
			}
			print PRIMER3 "\n";
		}
	}
}

close MFA;



