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
#  Module: Simulation::Simulate_NGS::random_insertion_set.pm
#  Purpose:
#  In:
#  Out:
#


###############################################
# Insertion simulation using multiple sequence
# sources like other genomes, transposable
# elements and tandem gene duplication
###############################################

package random_insertion_set;

our $VERSION = '1.0';

my $dbh;

sub new {
	my $self = shift;
	$self = {

	};
	bless $self;
	return $self;
}


###################################################################################################
# Creates random inserts for one chromosome using a reference chromosome, two other genomes, 
# transposon sequences and tandem gene duplication
sub create {

	my ($self, $ref_seq, $genome1_file, $genome2_file, $transposon_file, $annotation_file, $tandem_frequency, $total_bases, $outfile) = @_;
	my @nucleotides = ('A', 'T', 'C', 'G');
	my %ins = ();
	my %transposons = ();
	my $ref_len = length($ref_seq);


	### Read sequence from genome 1 (first source for insertion sequence)
	open( GENOME1, $genome1_file ) or die "Cannot open input file $genome1_file\n";
	my $genome_seq_1 = "";
	while( <GENOME1> ) {
		chomp;
		if($_ !~ />/) {
			$_ =~ s/N//g;
			$genome_seq_1 .= $_;
		}
	}


	### Read sequence from genome 2 (second source for insertion sequence)
	open GENOME2, $genome2_file or die "Cannot open input file $genome2_file\n";
	my $genome_seq_2 = "";
	while( <GENOME2> ) {
		chomp;
		if($_ !~ />/) {
			$_ =~ s/N//g;
			$genome_seq_2 .= $_;
		}
	}

	
	### Read fasta file of transposable elements (third source for insertion sequence
	open TE, $transposon_file or die "Cannot open input file $transposon_file\n";
	my $te_seq = "";
	my $te_count = -1;
	while( <TE> ) {
		chomp;
		if($_ =~ />/) {
			if($te_count != -1) {
				$transposons{$te_count} = $te_seq;
				$te_seq = "";
			}
			$te_count++;
		}
		else {
			$te_seq .= $_;
		}
	}


	### Create random set of insertions
	for( my $i = 0; $i < $total_bases; $i++ ) {
	
		# Randomize position and insertion sequence source
		my $pos = int(rand($ref_len)) + 1;
		my $select_sequence_source = int(rand(3));

		if( ! exists $ins{$pos} ) {

			# Randomize insertion length
			my $max_add    = 0;
			my $min_length = 1;
			my $length_type = int(rand(100));
				
			if($length_type < 1) {
				$max_add    = 50000;
				$min_length = 1000;
			}
			elsif($length_type < 4) {
				$max_add    = 900;
				$min_length = 100;
			}
			elsif($length_type < 10) {
				$max_add    = 95;
				$min_length = 5;
			}
			elsif($length_type < 30) {
				$max_add    = 3;
				$min_length = 2;
			}
			else {
				$max_add    = 0;
				$min_length = 1;
			}
			my $len = int(rand($max_add)) + $min_length;


			# Short insertions have random seq
			if($len <= 3) {
				my $ins_event = int(rand(4));
				$ins{$pos} = $nucleotides[$ins_event];

				if($len > 1) {
					my $ins_event = int(rand(4));
					$ins{$pos} .= $nucleotides[$ins_event];
					$i++;
				}

				if($len > 2) {
					my $ins_event = int(rand(4));
					$ins{$pos} .= $nucleotides[$ins_event];
					$i++;
				}
			}


			# Long insertions come from other sequence sources
			else {
				# Insert long seq from genome 1
				if($select_sequence_source == 0) {
					my $random_pos = int(rand(length($genome_seq_1)));
					$ins{$pos} = substr($genome_seq_1, $random_pos, $len);
					$i = $i + $len - 1;
				}

				# Insert long seq from genome 2
				elsif($select_sequence_source == 1) {
					my $random_pos = int(rand(length($genome_seq_2)));
					$ins{$pos} = substr($genome_seq_2, $random_pos, $len);
					$i = $i + $len - 1;
				}

				# Insertions originating from a transposable element. Complete transposon
				# sequence is inserted independent of randomized length
				elsif( $select_sequence_source == 2 ) {
				
					my $transposon_num = (int(rand($te_count)));
					$ins{$pos} = $transposons{$transposon_num};
					$i = $i + length($transposons{$transposon_num}) - 1;
				}
			}
		}
		else { $i--; }
	}


	### Tandemly duplicated genes
	if( $annotation_file ne "" && $tandem_frequency > 0 ) {
		my %gene = ();
		my $gene_count = 0;

		# Read gene annotation file (gff format)
		open GFF, $annotation_file or die "Cannot open input file $annotation_file\n";
		while( <GFF> ) {
			chomp;
			my ($tair_id, $beg, $end) = split("\t", $_);
			$gene{$beg} = $end;
			$gene_count++;
		}

		# Insert tandemly duplicated genes
		foreach my $beg ( sort {$a<=>$b} %gene ) {

			my $dice = int(rand(int($gene_count * $tandem_frequency)));
			if($dice == 1) {
				my $len = $gene{$beg} - $beg + 1;
				my $distance = int(rand(3000));
				my $tandem_beg = $gene{$beg} + $distance;
				my $seq = substr($ref_seq, $tandem_beg, $len);

				for(my $i=0; $i < $len; $i++) {
					if( (int(rand(100))) < 5 ) {
						my $mutated_base = &randomize_base( substr($seq, $i, 1) );
						substr( $seq, $i, 1, $mutated_base );
					}
				}
			}
		}
	}


	### Write insertions to file
	open OUTFILE, ">$outfile" or die "Cannot open output file $outfile\n";
	foreach my $ins_pos (sort {$a<=>$b} keys %ins) {
		print OUTFILE "$ins_pos\t" . $ins{$ins_pos} . "\n";
	}
}


####################################################
### Randomize new base
sub randomize_base {
	my ($self, $base) = @_;
	my @nucleotides = ('A', 'T', 'C', 'G');
	my $mutated_base = $base;
	while($mutated_base eq $base) {
		my $snp_event = int(rand(4));
		$mutated_base = $nucleotides[$snp_event];
	}
	return($mutated_base);
}


####################################################
## Load insertions from flat file
sub read {
	my ($self, $infile) = @_;
	my %insertion = ();

	open INFILE, $infile  or die "Cannot open input file $infile\n";
	while( <INFILE> ) {
		chomp;
		my ($pos, $seq)  = split(/\t/, $_);
		$insertion{$pos} = $seq;
	}

	return(%insertion);
}

1;

__END__

=head1 NAME

random_insertion_set

=head1 DESCRIPTION

Create random set of insertion for the genome of Arabidopsis thaliana

=head1 METHODS

=item C<$insertion_creator = new random_insertion_set()>;
=item C<%insertion = $insertion_creator->create($ref_file, $genome1_file, $genome2_file, $transposon_file, $annotation_file, $tandem_frequency, $total_bases, $outfile)>;
=item C<%insertion = $insertion_creator->read("$insertion_file")>;

=head1 GETTER/SETTER

=head1 AUTHOR

Stephan Ossowski <stephan.ossowski@tuebingen.mpg.de>

=head1 LICENCE AND COPYRIGHT

Copyright (C) 2005 by Max Planck Institute for Developmental Biology,
Tuebingen.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

# vim: ft=perl sw=4 ts=4 expandtab
