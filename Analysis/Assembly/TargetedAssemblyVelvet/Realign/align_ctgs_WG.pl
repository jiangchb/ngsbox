#! /usr/bin/perl
use strict;
use lib "$ENV{PGSP}/Assembly/Realign/";
use needle_alignment;
use DBI;

my $dbh;
my $database = "solexa";
connect_to_db();

my $ecotype = "Bur-0";
print STDERR "Selected ecotype = $ecotype\n";

open MFA, "> LCR_assembly_alignment.$ecotype.mfa";
open SNP, "> LCR_snp_table.$ecotype.txt";
open DEL, "> LCR_del_table.$ecotype.txt";
open INS, "> LCR_ins_table.$ecotype.txt";

my $max_contig = 1;

my $perc_recover = 0.80;
my $perc_foreign = 0.05;
my $perc_recover_flanking = 0.80;
my $min_reads_flanking = 10;
my $min_lo = 5;
my $min_align_length = 150;
my $trim_length = 6;
my $max_missing_nucs = 5;

my $debug_flag = 0;


my $usage = "perl align_reads_WG.pl contig_content.txt";

#system("sort -k4 $ARGV[0] > tmp");
#system("mv tmp $ARGV[0]");

open CTG, $ARGV[0] or die $usage;
open OUTPUT, ">".$ecotype.".WG_alignments.txt";

my $ctg_count = 0;

my $count_rejected_recover = 0;
my $count_rejected_flanking_recover = 0;
my $count_rejected_flanking_min_reads = 0;
my $count_rejected_foreign = 0;
my $count_rejected_no_mn_overlap = 0;
my $count_rejected_lo = 0;
my $count_rejected_align_length = 0;
my $count_rejected_missing_flanking = 0;

my $count_accepted = 0;
my $count_accepted_perfect = 0;
my $count_accepted_errors = 0;

my $count_accepted_seqnum = 0;
my $count_rejected_seqnum = 0;

while (my $line = <CTG>) {
	$ctg_count++;
	print $ctg_count, "\n" if $ctg_count%5000 == 0;

	my @a = split " ", $line;

	if ($a[2] ne "lo_contig") {

		my $ctg_id = $a[0];
		my $ctg_len = $a[1];
		my $ctg_seq = $a[15];

		my @b = split "-", $a[2];
		my $pcpr_chr = $b[0];
		my $pcpr_start = $b[1];
		my $pcpr_end = $b[2];

		my $pcpr_real_start = $a[4];
		my $pcpr_real_end = $a[5];
		my $pcpr_num_expected = $a[6];
		my $pcpr_num_recovered = $a[7];
		my $pcpr_num_foreign = $a[8] - $a[7];	
		my $pcpr_num_lo = $a[10];
	
		my $pcpr_num_expected_5prime = $a[11];
		my $pcpr_num_recovered_5prime = $a[12];
		my $pcpr_num_expected_3prime = $a[13];
		my $pcpr_num_recovered_3prime = $a[14];

		my $col_seq = "";

		if ($pcpr_real_start<$pcpr_start or $pcpr_real_end>$pcpr_end) {
			$count_rejected_seqnum++;
		}
		else {
			$count_accepted_seqnum++;

			#print STDERR $pcpr_num_recovered, "\t", $pcpr_num_expected, "\n";	
			if ($pcpr_num_recovered/$pcpr_num_expected >= $perc_recover) {
				if ($pcpr_num_expected_5prime >= $min_reads_flanking and $pcpr_num_expected_3prime >= $min_reads_flanking) {
					if ($pcpr_num_recovered_5prime/$pcpr_num_expected_5prime >= $perc_recover_flanking and $pcpr_num_recovered_3prime/$pcpr_num_expected_3prime >= $perc_recover_flanking) {
						if ($pcpr_num_foreign/($pcpr_num_foreign+$pcpr_num_recovered) <= $perc_foreign) {
							if ($pcpr_num_lo >= $min_lo) {

								my $q2 = "	SELECT  base  
										FROM 	seq_ref
										WHERE 	chromosome = $pcpr_chr and 
											position between $pcpr_start+$trim_length and $pcpr_end-$trim_length
										ORDER BY position";
					
								my $sth2 = $dbh->prepare($q2); $sth2->execute();
								while (my $res = $sth2->fetchrow_hashref()) {
									$col_seq .= $res->{base};
								}

								my $aligner = new needle_alignment();
								$aligner->needleman_wunsch($col_seq, substr($ctg_seq, $trim_length, length($ctg_seq)-(2*$trim_length)));
								$aligner->parse_fasta_alignment();

								my $alignment_length = $aligner->{end} - $aligner->{start} + 1;
								if ($alignment_length >= $min_align_length) {
									if ($aligner->{start} <= $max_missing_nucs and length($aligner->{align_seq1}) - $aligner->{end} <= $max_missing_nucs) {

										$count_accepted++;

										if ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and $aligner->{num_ins} == 0) {
											$count_accepted_perfect++;
										}
										else {
											$count_accepted_errors++;	
										}

										######### Screen Output #####################################
										if (1) {
											print "####################################################\n";
											print "Score:", $aligner->{similarity}, "\n";	
											print "PCPR:\t\t", $pcpr_chr, "\t", $pcpr_start+$trim_length, "\t", $pcpr_end-$trim_length, "\n";
											print "RealPCPR:\t", $pcpr_real_start, "\t", $pcpr_real_end, "\n";
											print "NumSNPS:", $aligner->{num_subs}, "\tNumDels:", $aligner->{num_del}, "\tNumIns:", $aligner->{num_ins}, "\n";
											print "NumRecovered:", $pcpr_num_recovered, "\tNumExpected:", $pcpr_num_expected, "\n";
											print "NumRecovered5Prime:", $pcpr_num_recovered_5prime, "\tNumExpected5Prime:", $pcpr_num_expected_5prime, "\n";
											print "NumRecovered3Prime:", $pcpr_num_recovered_3prime, "\tNumExpected3Prime:", $pcpr_num_expected_3prime, "\n";
											print "NumForeign:", $pcpr_num_foreign, "\tNumLO:", $pcpr_num_lo, "\n";
											print $aligner->{align_seq1}, "\t", length($aligner->{align_seq1}), "\n";
											print $aligner->{align_seq2}, "\t", length($aligner->{align_seq2}),"\n";
										}

										print OUTPUT $ecotype, "\t", $ctg_id, "\t", $pcpr_chr, "\t", $pcpr_start+$trim_length, "\t", $pcpr_end-$trim_length, "\n";

										######## MFA Out #############################################
										if (1) {
											print MFA ">LCR_", $count_accepted, " reference Col-0 ", $pcpr_chr, " ", $pcpr_start+$trim_length, " ", $pcpr_end-$trim_length, "\n";
											print MFA $aligner->{align_seq1}, "\n";
											print MFA ">LCR_", $count_accepted, " contig $ecotype \n";
                                                                                        print MFA $aligner->{align_seq2}, "\n";
											print MFA "\n";
										}

										####### Parse Polymorphismas with Col-Location ###############lignment_length >= $min_align_length
										
										parse_polymorphic_locations($ctg_id, $pcpr_chr, $pcpr_start+$trim_length, $pcpr_real_start, $pcpr_real_end, $aligner->{align_seq1}, $aligner->{align_seq2});
										$debug_flag = 1;

									}
									else {
										$count_rejected_missing_flanking++;
									}
								}
								else {
									$count_rejected_align_length++;
								}
							}
							else {
								$count_rejected_lo++;
							}
						}
						else {
							$count_rejected_foreign++;
						}
					}
					else {
						$count_rejected_flanking_recover++;
					}
				}
				else {
					$count_rejected_flanking_min_reads++;
				}
			}
			else {
				$count_rejected_recover++;
			}
		}
	}

	if ($debug_flag == 2) {
		close MFA;
		close SNP;
		close DEL;
		close INS;
		exit(1);
	}

}

close MFA;
close SNP;
close DEL;
close INS;
close OUTPUT;

if (1) {
        print "Accepted:numseq:", $count_accepted_seqnum, "\n";
        print "Rejected:numseq:", $count_rejected_seqnum, "\n";
        print "\n";
	print "Accepted:", $count_accepted, "\n";
	print "Accepted:perfect:", $count_accepted_perfect, "\n";
	print "Accepted:errors:", $count_accepted_errors, "\n";
	print "\n";
	print "Rejected:foreign:", $count_rejected_foreign, "\n";
	print "Rejected:recover:", $count_rejected_recover, "\n";
	print "Rejected:flankingrecover:", $count_rejected_flanking_recover, "\n";
	print "Rejected:flankingnumreads:", $count_rejected_flanking_min_reads, "\n";
	print "Rejected:nooverlap:", $count_rejected_no_mn_overlap, "\n"; 
	print "Rejected:lo:", $count_rejected_lo, "\n";
	print "rejected:alignmentlength:", $count_rejected_align_length, "\n";
	print "rejected:missingnucs:", $count_rejected_missing_flanking, "\n";
}


sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

sub parse_polymorphic_locations {
	my ($count_ctg, $chr, $start, $pcpr_real_start, $pcpr_real_end, $seq1, $seq2) = @_;

	my $pos_chr = $start;
	my $pos_str = 0;
	my $flag = 0;

	while (substr($seq1, $pos_str, 1) eq "-") {
		if (substr($seq2, $pos_str, 1) ne "-" and substr($seq2, $pos_str, 1) ne "N") {
			$flag = 1;
		}
		$pos_str++;
	}

	while (substr($seq2, $pos_str, 1) eq "-" or substr($seq2, $pos_str, 1) eq "N") {
		if (substr($seq1, $pos_str, 1) ne "-") {
			$pos_chr++;
		}
		$pos_str++;
	}

	my $insertion_flag = 0;
	my $insertion_start;
	my $insertion_seq = "";

	my $deletion_flag = 0;
	my $deletion_start;
	my $deletion_seq = "";

	for (; $pos_str < length($seq1); $pos_str++) {

		my $base1 = substr($seq1, $pos_str, 1);
		my $base2 = substr($seq2, $pos_str, 1);

		# End of insertion
		if ($insertion_flag == 1 and $base1 ne "-") {
			if (($pos_chr >= $pcpr_real_start and $pos_chr <= $pcpr_real_end) or ($insertion_start >= $pcpr_real_start and $insertion_start <= $pcpr_real_end)) {
				print INS $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $insertion_start, "\t", $pos_chr, "\t", "LCR", "\t", length($insertion_seq), "\t", $insertion_seq, "\n";
			}
			else {
				print INS $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $insertion_start, "\t", $pos_chr, "\t", "FLANKING", "\t", length($insertion_seq), "\t", $insertion_seq, "\n";
			}
			$insertion_flag = 0;
			$insertion_seq = "";
		}

		# End of deletion
		if ($deletion_flag == 1 and $base2 ne "-") {
                        if (($pos_chr >= $pcpr_real_start and $pos_chr <= $pcpr_real_end) or ($deletion_start >= $pcpr_real_start and $deletion_start <= $pcpr_real_end)) {
                                print DEL $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $deletion_start, "\t", $pos_chr, "\t", "LCR", "\t", length($deletion_seq), "\t", $deletion_seq, "\n";
                        }
                        else {
                                print DEL $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $deletion_start, "\t", $pos_chr, "\t", "FLANKING", "\t", length($deletion_seq), "\t", $deletion_seq, "\n";
                        }
                        $deletion_flag = 0;
                        $deletion_seq = "";
                }		


		if ($base1 eq "-") {
			if ($insertion_flag != 1) {
				$insertion_flag = 1;
				$insertion_start = $pos_chr;
			}
			$insertion_seq .= $base2;
		}
		elsif ($base2 eq "-") {
			if ($deletion_flag != 1) {
				$deletion_flag = 1;
				$deletion_start = $pos_chr;
			}
			$deletion_seq .= $base1;
			$pos_chr++;
		}
		elsif ($base1 ne $base2 and $base1 ne "-" and $base2 ne "-") {
			if ($pos_chr >= $pcpr_real_start and $pos_chr <= $pcpr_real_end) {
				print SNP $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $pos_chr, "\t", $base1, "\t", $base2, "\t", "LCR", "\n";
			}
			else {
				print SNP $ecotype, "\t", $count_ctg, "\t", $chr, "\t", $pos_chr, "\t", $base1, "\t", $base2, "\t", "FLANKING", "\n";
			}
			$pos_chr++;
		}
		else {
			$pos_chr++;
		}
		
	
	}


}

