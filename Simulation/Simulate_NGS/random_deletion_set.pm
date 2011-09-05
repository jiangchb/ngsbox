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
#  Module: Simulation::Simulate_NGS::random_deletion_set.pm
#  Purpose:
#  In:
#  Out:
#


##############################################
# Deletion simulation
##############################################

package random_deletion_set;

our $VERSION = '0.9';


sub new {
	my $self = shift;
	$self = {

	};
	bless $self;
	return $self;
}

### Creates a simulated deletion set
sub create {

	my ($self, $ref_len, $total_bases, $outfile) = @_;
	my %deletion = ();
	my %del_conc = ();

	### Create random set of deletion
	my $i = 0;
	while( $i < $total_bases ) {
	
		# Randomize chromosome and position
		my $beg = int(rand($ref_len)) + 1;
	
		# Randomize length
		my $max_add     = 0;
		my $min_length  = 1;
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
			$max_add    = 99;
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

		my $del_length = int(rand($max_add)) + $min_length;

		# Check if overlapping with previous deletion
		for(my $j = $beg; $j < ($beg + $del_length); $j++) {
			if(! exists $deletion{$j}) {
				$deletion{$j} = 1;
				$i++;
			}
		}
	}


	### Write deletions to file (concatenate adjacent gaps)
	my $beg = -10;
	my $end  = -10;
	open OUTFILE, ">$outfile"  or die "Cannot open output file $outfile\n";
	foreach my $pos (sort {$a<=>$b} keys %deletion) {

		if( ($pos != $end + 1) ) {
			if($beg != -10) {
				print OUTFILE "$beg\t$end\n";
				$del_conc{$beg} = $end;
			}
			$beg = $pos;
		}
		$end = $pos;
	}
	print OUTFILE "$beg\t$end\n";
	$del_conc{$beg} = $end;

	### Return hash of deletions
	return(%del_conc);
}


### Read a set of deletions
sub read {
	my ($self, $infile) = @_;
	open INFILE, "$infile"  or die "Cannot open input file $infile\n";
	my %deletion = ();

	while( <INFILE> ) {
		chomp($_);
		my ($beg, $end) = split("\t", $_);
		$deletion{$beg} = $end;
	}

	return(%deletion);
}

1;

__END__

=head1 NAME

random_deletion_set

=head1 DESCRIPTION

Create random set of deletion for the genome of Arabidopsis thaliana

=head1 METHODS

=item C<$deletion_creator = new random_deletion_set()>;
=item C<%deletion = $deletion_creator->create("$deletion_file")>;
=item C<%dels = $deletion_creator->read("$deletion_file")>;

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
