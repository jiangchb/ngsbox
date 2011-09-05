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
#  Module: Simulation::Simulate_NGS::random_polymorphisms.pl
#  Purpose:
#  In:
#  Out:
#


our $VERSION = '1.0';


use Getopt::Long;
use FindBin;
use lib $FindBin::Bin;

use random_deletion_set;
use random_insertion_set;
use random_snp_set;


### User parameters
my $reference_file   = "";
my $genome1_file     = "";
my $genome2_file     = "";
my $transposon_file  = "";
my $annotation_file  = "";
my $deletion_file    = "";
my $insertion_file   = "";
my $snp_file         = "";
my $snp_num          = 600000;
my $del_num          = 1000000;
my $ins_num          = 1000000;
my $tandem_frequency = 0.0;

my %CMD;
GetCom();


### Read reference chromosome
open REF, $reference_file or die "Cannot open input file $reference_file\n";
my $ref_seq = "";
while( <REF> ) {
	chomp;
	if($_ !~ /[N>]/) {
		$ref_seq .= $_;
	}
}
my $ref_len = length($ref_seq);


# Create random set of deletions
my $deletion_creator = new random_deletion_set();
my %deletion = $deletion_creator->create($ref_len, $del_num, $deletion_file);


# Create random set of insertions
my $insertion_creator = new random_insertion_set();
$insertion_creator->create($ref_seq, $genome1_file, $genome2_file, $transposon_file, 
			$annotation_file, $tandem_frequency, $ins_num, $insertion_file);


# Create random set of SNPs
my $snp_creator = new random_snp_set();
$snp_creator->create($ref_seq, $snp_num, $snp_file, %deletion);


exit(0);


### Read command line parameters
sub GetCom {

	my @usage = ("$0

Mandatory:
--reference    STRING    Reference chromosome
--genome1      STRING    Chromosome from different species in fasta format
--genome2      STRING    Chromosome from different species in fasta format
--transposon   STRING    Transposons in fasta format
--delfile      STRING    Outfile for deletions
--insfile      STRING    Outfile for insertions
--snpfile      STRING    Outfile for SNPs

Optional:
--snp          INT       Number of SNP to simulated                  (default: 600,000)
--del          INT       Total number of nucleotides to be deleted   (default: 1MB)
--ins          INT       Total number of nucleotides to be inserted  (default: 1MB)
--annotation   STRING    Gene annotation in GFF format               (default: none)
--frequency    FLOAT     Frequency of tandemly duplicated genes      (default: 0.0)

\n");

	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "reference=s","genome1=s", "genome2=s", "transposon=s", "delfile=s", "insfile=s", 
			"snpfile=s", "snp=s","del=s","ins=s","annotation=s","frequency=s");

	die("Please specify reference infile\n") unless defined($CMD{reference});
	die("Please specify genome1 infile\n") unless defined($CMD{genome1});
	die("Please specify genome2 infile\n") unless defined($CMD{genome2});
	die("Please specify transposon infile\n") unless defined($CMD{transposon});
	die("Please specify deletion outfile\n") unless defined($CMD{delfile});
	die("Please specify insertion outfile\n") unless defined($CMD{insfile});
	die("Please specify SNP outfile\n") unless defined($CMD{snpfile});


	$reference_file  = $CMD{reference};
	$genome1_file    = $CMD{genome1};
	$genome2_file    = $CMD{genome2};
	$transposon_file = $CMD{transposon};
	$deletion_file   = $CMD{delfile};
	$insertion_file  = $CMD{insfile};
	$snp_file        = $CMD{snpfile};

	if (defined($CMD{snp}))        { $snp_num = $CMD{snp}; }
	if (defined($CMD{del}))        { $del_num = $CMD{del}; }
	if (defined($CMD{ins}))        { $ins_num = $CMD{ins}; }
	if (defined($CMD{annotation})) { $annotation_file = $CMD{annotation}; }
	if (defined($CMD{frequency}))  { $tandem_frequency = $CMD{frequency}; }
}

__END__

=head1 NAME

create_random_polymorphisms

=head1 DESCRIPTION

randomly create SNPs, deletions and insertions for sequencing simulations

=head1 METHODS

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
#
