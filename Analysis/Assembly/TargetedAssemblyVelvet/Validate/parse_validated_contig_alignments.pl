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
#  Module: Analysis::Assembly::TargetedAssemblyVelvet::Validate::parse_validated_contig_alignments.pl
#  Purpose:
#  In:
#  Out:
#

use lib "$ENV{PGSP}/Assembly/Realign/";
use needle_alignment;

open CONTIG_IDS, $ARGV[0] or die "no file\n";
open ALIGN_BUR_COL, $ARGV[1] or die "no file";

my %ID = ();

my $id = "";
my $seq = "";

while (my $line = <CONTIG_IDS>) {
	chomp($line);
	$ID{$line} = $line;	
}


my $count_del_3 = 0;
my $del_max = 0;
my $count_ins_3 = 0;
my $ins_max = 0;

while (my $line = <ALIGN_BUR_COL>) {
	chomp($line);
        if ($line =~ /reference/) {
		my @a = split " ", $line;
		my $colseq = <ALIGN_BUR_COL>;
                my $junk = <ALIGN_BUR_COL>;
                my $burseq = <ALIGN_BUR_COL>;
                my $junk2 = <ALIGN_BUR_COL>;
		if (defined($ID{substr($a[0], 5, length($a[0])-5)})) {
			$colseq =~ s/-//g;
			$burseq =~ s/-//g;
			my $aligner = new needle_alignment();
		        $aligner->needleman_wunsch($colseq, $burseq);
		        $aligner->parse_fasta_alignment();
	
			print "NumSNPS:", $aligner->{num_subs}, "\tNumDels:", $aligner->{num_del}, "\tNumIns:", $aligner->{num_ins}, "\n";
			foreach my $key (keys %{$aligner->{del}}) {
				if (length($aligner->{del}{$key}) > 3) {
					$count_del_3++;
				}
				if (length($aligner->{del}{$key}) > $del_max) {
					$del_max = length($aligner->{del}{$key});
				}
			}
			foreach my $key (keys %{$aligner->{ins}}) {
				if (length($aligner->{ins}{$key}) > 3) {
                                        $count_ins_3++;
                                }
                                if (length($aligner->{ins}{$key}) > $ins_max) {
                                        $ins_max = length($aligner->{ins}{$key});
                                }
			}
		}
		
	}
}

print "Num deletions > 3: ", $count_del_3, "\n";
print "Max length deletions: ", $del_max, "\n";
print "Num insertions > 3:", $count_ins_3, "\n";
print "Max length insertion: ", $ins_max, "\n";


# 
# sub align {
# 	my ($id, $seq1, $seq2) = @_;
# 
# 
# 	my $aligner = new needle_alignment();
# 	$aligner->needleman_wunsch($seq1, $seq2);
# 	$aligner->parse_fasta_alignment();
# 
# 	#my $align = needle_alignment::needleman_wunsch($seq1, $seq2);
# 	#my ($snp_ref, $del_ref, $ins_ref) = needle_alignment::parse_alignment($align, 0);
# 
# 	if (not ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0)) {
# 		print "\n\n";
# 		print ">$id\n";
# 		print $aligner->{align_seq1}, "\n";
#         	print $aligner->{align_seq2}, "\n";
# 		for (my $i = 0; $i<length($aligner->{align_seq1}); $i++) {
# 			if (uc(substr($aligner->{align_seq2}, $i, 1)) ne uc(substr($aligner->{align_seq1}, $i, 1))) {
# 				print "X";
# 			}
# 			else {
# 				print " ";
# 			}
# 		}
# 		print "\n";
# 
# 		print "NumSNPS:", $aligner->{num_subs}, "\tNumDels:", $aligner->{num_del}, "\tNumIns:", $aligner->{num_ins}, "\n";
# 	}
# 
# 
# 	$count_seq++;		
# 	if ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0) {
# 		$count_perfects++;
# 	}
# 	elsif ($aligner->{num_subs} < 3 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0) {
# 		$count_nearly_perfect++;
# 	}
# 
# }


