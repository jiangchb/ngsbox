#! /usr/bin/perl

my $usage = "perl parse_reads_in_ctgs.pl ctg2reads.txt pcpr2reads.txt";
open CTG, $ARGV[0] or die $usage;
open PCPR, $ARGV[1] or die $usage;

my %READ2PCPR = ();
my %READ2PCPR_5prime = ();
my %READ2PCPR_3prime = ();
my %PCPR2NUM = ();
my %PCPR2NUM_5prime = ();
my %PCPR2NUM_3prime = ();
my %PCPR2REAL_START = ();
my %PCPR2REAL_END = ();
while (<PCPR>) {
	my @a = split " ", $_;;

	my @r = split ",", $a[5];
	for (my $i = 0; $i < @r; $i++) {
		if (not defined($READ2PCPR{$r[$i]})) {
			$READ2PCPR{$r[$i]} = $a[0]."-".$a[1]."-".$a[2];
		}
		else {
			$READ2PCPR{$r[$i]} = "m"; # flag reads with the same id -> useless;
		}
	}
	
	my @r_5prime = split ",", $a[6];
        for (my $i = 0; $i < @r_5prime; $i++) {
		if ($READ2PCPR{$r_5prime[$i]} ne "m") {
	        	$READ2PCPR_5prime{$r_5prime[$i]} = $a[0]."-".$a[1]."-".$a[2];
		}
		else {
			$READ2PCPR_5prime{$r_5prime[$i]} = "m";
		}
        }

	my @r_3prime = split ",", $a[7];
        for (my $i = 0; $i < @r_3prime; $i++) {
		if ($READ2PCPR{$r_3prime[$i]} ne "m") {
	                $READ2PCPR_3prime{$r_3prime[$i]} = $a[0]."-".$a[1]."-".$a[2];
		}
		else {
			$READ2PCPR_3prime{$r_3prime[$i]} = "m";
		}
        }

	$PCPR2NUM{$a[0]."-".$a[1]."-".$a[2]} = @r + 0;
	$PCPR2NUM_5prime{$a[0]."-".$a[1]."-".$a[2]} = @r_5prime + 0;
	$PCPR2NUM_3prime{$a[0]."-".$a[1]."-".$a[2]} = @r_3prime + 0;
	$PCPR2REAL_START{$a[0]."-".$a[1]."-".$a[2]} = $a[3];
	$PCPR2REAL_END{$a[0]."-".$a[1]."-".$a[2]} = $a[4];
}



while (<CTG>) {
	my @a = split " ";

	# go through reads of the ctg
	my $ctg_id = $a[0];
	my $ctg_len = $a[1];
	my $ctg_seq = $a[2];

	my $count_lo = 0;
	my $count_unvalid = 0;
	my $count_all_pcpr_reads = 0;
	my %pcpr = ();

	for (my $i = 3; $i < @a; $i++) {
		if ($a[$i] =~ /-run_pc_pr/) {
			$a[$i] =~ s/-run_pc_pr/-pc/g;
			my $pcpr_id = $READ2PCPR{$a[$i]};
			if ($pcpr_id ne "m") {
				$pcpr{$pcpr_id}++;
				$count_all_pcpr_reads++;
			}
			else {
				$count_unvalid++;
			}
		}
		else {
			$count_lo++;
		}
	}

	# Set pcpr_id with the most reads in the actual contig as the target pcpr_id
	# else it is a left over contig
	my $recovered_read_num = 0;
	my $pcpr_id = "lo_contig";
	my $pcpr_len = 0;
	my $max = 0;
	foreach my $key (keys %pcpr) {
		if ($max < $pcpr{$key}) {
			$recovered_read_num = $pcpr{$key};
			$pcpr_id = $key;
			my @a = split "-", $pcpr_id;
			$pcpr_len = $a[2] - $a[1] + 1;	
		}
	}	

	# Set variable specific to the detected pcpr
	my $exp_recovered_read_num = 0;
	my $exp_recovered_read_num_5prime = 0;
	my $exp_recovered_read_num_3prime = 0;
	my $recovered_read_num_5prime = 0;
	my $recovered_read_num_3prime = 0;
	my $pcpr_real_start = -1;
	my $pcpr_real_end = -1;

	if ($pcpr_id ne "lo_contig") {
		$exp_recovered_read_num = $PCPR2NUM{$pcpr_id};
		$exp_recovered_read_num_5prime = $PCPR2NUM_5prime{$pcpr_id};
		$exp_recovered_read_num_3prime = $PCPR2NUM_3prime{$pcpr_id};
		$pcpr_real_start = $PCPR2REAL_START{$pcpr_id};
		$pcpr_real_end = $PCPR2REAL_END{$pcpr_id};
		for (my $i = 3; $i < @a; $i++) {
			if ($a[$i] =~ /-pc/) {
				if (defined($READ2PCPR_5prime{$a[$i]}) and $READ2PCPR_5prime{$a[$i]} eq $pcpr_id) {
					$recovered_read_num_5prime++;
				}
				if (defined($READ2PCPR_3prime{$a[$i]}) and $READ2PCPR_3prime{$a[$i]} eq $pcpr_id) {
                                        $recovered_read_num_3prime++;
                                }
			}	
		} 
	}	

	print $ctg_id, "\t", $ctg_len, "\t", $pcpr_id, "\t", $pcpr_len, "\t", $pcpr_real_start, "\t", $pcpr_real_end, "\t", $exp_recovered_read_num, "\t", $recovered_read_num, "\t", $count_all_pcpr_reads, "\t", $count_unvalid, "\t", $count_lo, "\t", $exp_recovered_read_num_5prime, "\t", $recovered_read_num_5prime, "\t", $exp_recovered_read_num_3prime, "\t", $recovered_read_num_3prime, "\t", $ctg_seq, "\n";
}

