#!/usr/bin/perl 

##############################################
# Simulates or reads a set of SNPs including
# dense SNP cluster
##############################################

package random_snp_set;

our $VERSION = '0.09';

use strict;
use DBI;

my $dbh = "";

sub new {
	my $self = shift;
	$self = {

	};
	bless $self;
	return $self;
}


### Creates a simulated SNP set
sub create {

	my ($self, $ref_seq, $snp_num, $outfile, %del) = @_;
	my $chr_len = length($ref_seq);
	my %snps = ();
	my %deletion = ();

	### Read deletion by single position
	foreach my $beg ( sort {$a<=>$b} keys %del ) {
		for(my $i = $beg; $i <= $del{$beg}; $i++) {
			$deletion{$i} = 1;
		}
	}
	
	open OUTFILE, ">$outfile" or die "Cannot open output file $outfile\n";

	# Simulate SNPs
	for( my $i = 0; $i < $snp_num; $i++ ) {
	
		# Randomize chromosome and position
		my $pos = int(rand($chr_len)) + 1;

		# Write SNP
		if( (! exists $snps{$pos}) && (! exists $deletion{$pos}) ) {
			my $ref_base = substr($ref_seq, $pos - 1, 1);
			$snps{$pos} = &randomize_base($ref_base);

			# Highly polymorphic region
			my $hpr = int(rand(5));
			if($hpr == 2) {
				my $close_snp = int(rand(60));
				my $close_pos = $pos - 30 + $close_snp;
				if(! exists $snps{$close_pos}) {
					my $ref_base = substr($ref_seq, $close_pos - 1, 1);
					$snps{$close_pos} = &randomize_base($close_pos);
					$i++;
				}
			
				while( int(rand(3)) == 1 ) {
					$close_snp = int(rand(60));
					$close_pos = $pos - 30 + $close_snp;
					if(! exists $snps{$close_pos}) {
						my $ref_base = substr($ref_seq, $close_pos - 1, 1);
						$snps{$close_pos} = &randomize_base($close_pos);
						$i++;
					}
				}
			}
		}
		else { $i--; }
	}

	# Write SNPs to file
	foreach my $snp_pos (sort {$a<=>$b} keys %snps) {
		print OUTFILE "$snp_pos\t" . $snps{$snp_pos} . "\n";
	}
}


### Returnns random nucleotide different from the reference base
sub randomize_base {
	my $ref_base = shift;
	my @nucleotides = ('A', 'T', 'C', 'G');

	# Randomize new base
	my $snp_base = $ref_base;
	while($snp_base eq $ref_base) {
		my $snp_event = int(rand(4));
		$snp_base = @nucleotides[$snp_event];
	}
	return($snp_base);
}



sub read {
	my ($self, $infile) = @_;
	open INFILE, $infile  or die "Cannot open input file $infile\n";
	my %snps = ();

	while( <INFILE> ) {
		chomp($_);
		my ($pos, $snp_base) = split(/\t/, $_);
		$snps{$pos} = $snp_base;
	}
	return(%snps);
}

1;
__END__

=head1 NAME

random_snp_set

=head1 DESCRIPTION

Simulates a set of SNPs for one chromosome

=head1 METHODS

=item C<$snp_creator = new random_snp_set()>;
=item C<%snps = $snp_creator->create("$snp_file")>;
=item C<%snps = $snp_creator->read("$snp_file")>;

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
