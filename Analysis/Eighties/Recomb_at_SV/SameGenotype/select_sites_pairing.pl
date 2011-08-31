#! /usr/bin/perl
use strict;

my $min_freq = 5;
my $max_freq = 75;

my $usage = "$0 SV_SNP_pattern.txt\n";
open FILE, shift or die $usage;

my $count = 0;
while (my $pat_line = <FILE>) {
	$count++;
	#print STDERR $count, "\n" if $count%100 == 0;

	my $sv_line = <FILE>;
	my $snp_line = <FILE>;
	my $fuzzy_snp_line = <FILE>;

	chomp($pat_line);
	chomp($sv_line);
	chomp($snp_line);
	chomp($fuzzy_snp_line);
	
	#print STDERR "#############\n";
	#print STDERR ">", $pat_line;
	#print STDERR ">", $sv_line;
	#print STDERR ">", $snp_line;
	
	my @a = split " ", $pat_line;
	my @c = split ",", $a[1];

	if ($a[0] >= $min_freq and $a[0] <= $max_freq) {
		my @svs = split "#", $sv_line;
		my @snps = split "#", $snp_line;
		my @fuzzy_snps = split "#", $fuzzy_snp_line;
		SV: foreach my $sv (@svs) {
			if ($sv ne "") {

				my @b = split "-", $sv;
				my $sv_chr = int($b[0] / 100000000);
				my $sv_begin = $b[0]%100000000;
				my $sv_end = $b[1]%100000000;

				my @partner = ();
				my @partner_fuzzy = ();

				foreach my $snp (@snps) {
					if ($snp ne "") {
						my $snp_chr = int($snp / 100000000);
						my $snp_pos = $snp%100000000;
						if ($snp_chr != $sv_chr) {
							push @partner,  $sv_chr."\t".$sv_begin."\t".$sv_end."\t".$snp_chr."\t".$snp_pos."\t".scalar(@c)."\t".$a[1]."\n";
						}
					}
				}

				if (@partner+0 == 0 and 1==0) {
					foreach my $snp (@fuzzy_snps) {
        	                                if ($snp ne "") {
                	                                my $snp_chr = int($snp / 100000000);
                        	                        my $snp_pos = $snp%100000000;
                                	                if ($snp_chr != $sv_chr) {
								push @partner_fuzzy, $sv_chr."\t".$sv_begin."\t".$sv_end."\t".$snp_chr."\t".$snp_pos."\t".scalar(@c)."\t".$a[1]."\n";
							}
                                                }
                                        }
                                }

				my $res;
				if (@partner+0 != 0) {
					$res = \@partner;
				}
				elsif (@partner_fuzzy+0 != 0) {
					$res = \@partner_fuzzy;
				}
				else {
					next SV;
				}

				my $r = int(rand(@{$res}+0));
				print ${$res}[$r];

			}
		}
	}
}



