#! /usr/bin/perl
use strict;
use lib "$ENV{PGSP}/Assembly/Realign/";
use needle_alignment;
use DBI;

my $dbh;
my $database = "solexa";
connect_to_db();

my $ecotype = "Bur-0";

print "Selected ecotype = '$ecotype'\n";

my $max_contig = 1;

my $perc_recover = 0.80;
my $perc_foreign = 0.05;
my $perc_recover_flanking = 0.80;
my $min_reads_flanking = 10;
my $min_lo = 5;
my $min_align_length = 100;
my $trim_length = 6;

my $max_missing_nucs = 5;

my $usage = "perl align_reads_MN.pl contig_content.txt";

#system("sort -k4 $ARGV[0] > tmp");
#system("mv tmp $ARGV[0]");

open CTG, $ARGV[0] or die $usage;
open OUTPUT, ">".$ecotype.".MN_alignments.txt";

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
my $count_accepted_few_snps = 0;

my $count_accepted_seqnum = 0;
my $count_rejected_seqnum = 0;



while (my $line = <CTG>) {
	$ctg_count++;
	print $ctg_count, "\n" if $ctg_count%1000 == 0;

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

		my $frag_id = "";
		my $frag_seq = "";
		my $frag_start;
		my $frag_end;
		my $target_seq = "";
		my $target_seq_start;
		my $target_seq_end;	

		if ($pcpr_real_start<$pcpr_start or $pcpr_real_end>$pcpr_end) {
			$count_rejected_seqnum++;
		}
		else {
			$count_accepted_seqnum++;

			my $q2 = "	SELECT position, end, $pcpr_start-position as target_start, $pcpr_end-position as target_end, sequence 
					FROM mn_fragments a, mn_sequence b 
					WHERE 	a.ecotype = '$ecotype' and b.ecotype = '$ecotype' and 
						a.frag_id = b.frag_id and chromosome = $pcpr_chr and 
						position+num_ns_begin < $pcpr_real_start and end-num_ns_end > $pcpr_real_end";

			my $sth2 = $dbh->prepare($q2); $sth2->execute();
			while (my $res = $sth2->fetchrow_hashref()) {
				$frag_id = $res->{frag_id};
				$frag_seq = $res->{sequence};
				$frag_start = $res->{position};
				$frag_end = $res->{end};

				#$target_seq_start = $res->{target_start};
				#$target_seq_end = $res->{target_end};
				#$target_seq = substr($frag_seq, $target_seq_start, $target_seq_end-$target_seq_start+1);
			}


			if ($frag_seq ne "") {
				if ($pcpr_num_recovered/$pcpr_num_expected >= $perc_recover) {
					if ($pcpr_num_expected_5prime >= $min_reads_flanking and $pcpr_num_expected_3prime >= $min_reads_flanking) {
						if ($pcpr_num_recovered_5prime/$pcpr_num_expected_5prime >= $perc_recover_flanking and $pcpr_num_recovered_3prime/$pcpr_num_expected_3prime >= $perc_recover_flanking) {
							if ($pcpr_num_foreign/($pcpr_num_foreign+$pcpr_num_recovered) <= $perc_foreign) {
								if ($pcpr_num_lo >= $min_lo) {

									my $aligner = new needle_alignment();
						
									print "#################################################### $ctg_id $pcpr_num_lo\n";

									$aligner->needleman_wunsch($frag_seq, substr($ctg_seq, $trim_length, length($ctg_seq)-(2*$trim_length)));
									$aligner->parse_fasta_alignment();
				
		                                        	        #print $aligner->{align_seq1}, "\n";
        		                                        	#print $aligner->{align_seq2}, "\n";

									my $alignment_length = $aligner->{end} - $aligner->{start} + 1;
									if ($alignment_length >= $min_align_length) {
										print STDERR "##################################\n";
										print STDERR $frag_start, " ", $frag_end, "\n";
										print STDERR $pcpr_start, " ", $pcpr_end, "\n";
										print STDERR $aligner->{start}, " ", $aligner->{end}, "\n";
										#if 	(
										#	($frag_start > $pcpr_start or ($frag_start <= $pcpr_start and $aligner->{start} <= $max_missing_nucs)) 
										#	and 
										#	($frag_end < $pcpr_end or ($frag_end >= $pcpr_end and length($aligner->{align_seq1}) - $aligner->{end} <= $max_missing_nucs))
										# 	)
										#{

											$count_accepted++;

											#print "Number of SNPS:>".$aligner->{num_subs}."<\n";
       			                        		                        #print "Number of Deletions:>".$aligner->{num_del}."<\n";
       				                        	                        #print "Number of Insertions:>".$aligner->{num_ins}."<\n";

											if ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and $aligner->{num_ins} == 0) {
												$count_accepted_perfect++;
											}
											elsif ($aligner->{num_del} == 0 and $aligner->{num_ins} == 0 and $aligner->{num_subs} < 2) {
												$count_accepted_few_snps++;	
											}
	
											#if (not(($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and $aligner->{num_ins} == 0))) {
												print "####################################################\n";
												print "Score:", $aligner->{similarity}, "\n";
												print "Fragment:\t", $frag_start, "\t", $frag_end, "\n";	
												print "PCPR:\t\t", $pcpr_start, "\t", $pcpr_end, "\n";
												print "RealPCPR:\t", $pcpr_real_start, "\t", $pcpr_real_end, "\n";
												print "NumSNPS:", $aligner->{num_subs}, "\tNumDels:", $aligner->{num_del}, "\tNumIns:", $aligner->{num_ins}, "\n";
												print "NumRecovered:", $pcpr_num_recovered, "\tNumExpected:", $pcpr_num_expected, "\n";
												print "NumRecovered5Prime:", $pcpr_num_recovered_5prime, "\tNumExpected5Prime:", $pcpr_num_expected_5prime, "\n";
												print "NumRecovered3Prime:", $pcpr_num_recovered_3prime, "\tNumExpected3Prime:", $pcpr_num_expected_3prime, "\n";
												print "NumForeign:", $pcpr_num_foreign, "\tNumLO:", $pcpr_num_lo, "\n";
												print $aligner->{align_seq1}, "\t", length($aligner->{align_seq1}), "\n";
				                                               	                print $aligner->{align_seq2}, "\t", length($aligner->{align_seq2}),"\n";
												print $frag_seq, "\n", $ctg_seq, "\n";

												print OUTPUT $ecotype, "\t", $ctg_id, "\t", $pcpr_chr, "\t", $pcpr_start+$trim_length, "\t", $pcpr_end-$trim_length, "\t", $pcpr_chr ,"\t", $frag_start, "\t", $frag_end, "\n";
											#}
										#}
										#else {
										#	$count_rejected_missing_flanking++;
										#}
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
			else {
				$count_rejected_no_mn_overlap++;
			}
		}
	}
}

if (1) {
        print "Accepted:numseq:", $count_accepted_seqnum, "\n";
        print "Rejected:numseq:", $count_rejected_seqnum, "\n";
        print "\n";
	print "Accepted:", $count_accepted, "\n";
	print "Accepted:perfect:", $count_accepted_perfect, "\n";
	print "Accepted:few_snps:", $count_accepted_few_snps, "\n";
	print "\n";
	print "Rejected:foreign:", $count_rejected_foreign, "\n";
	print "Rejected:recover:", $count_rejected_recover, "\n";
	print "Rejected:flankingrecover:", $count_rejected_flanking_recover, "\n";
	print "Rejected:flankingnumreads:", $count_rejected_flanking_min_reads, "\n";
	print "Rejected:nooverlap:", $count_rejected_no_mn_overlap, "\n"; 
	print "Rejected:lo:", $count_rejected_lo, "\n";
	print "rejected:alignmentlength:", $count_rejected_align_length, "\n";
	print "rejected:missingflanking:", $count_rejected_missing_flanking, "\n";
}

close OUTPUT;

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



