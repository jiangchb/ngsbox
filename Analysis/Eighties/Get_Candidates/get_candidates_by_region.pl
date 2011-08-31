#! /usr/bin/perl
use strict;
use DBI;

my $usage = "$0 quality_reference.txt quality_variants.txt chr begin end\n";

my $file = shift or die $usage;
my $var = shift or die $usage;
my $chr = shift or die $usage;
my $begin = shift or die $usage;
my $end = shift or die $usage;

my %SNP = ();
my %REF = ();

open OUT, ">multiple_alignment.txt";

########################################################################

open FILE, $file or die "Cannot open file: ".$file."\n";

while (<FILE>) {
	my @a = split " ";
	if ($chr == $a[1] and $a[2] >= $begin and $a[2] <= $end) {
		if ($a[5] >= 25) {
			$REF{$a[2]} = $a[4];
		}
	}
	elsif ($chr < $a[1] or ($chr == $a[1] and $a[2] > $end)) {
		last;
	}
}

close FILE;

print STDERR "got ref: ", scalar(keys %REF), "\n";

########################################################################

open FILE, $var or die "Cannot open file: ".$var."\n";

while (<FILE>) {
        my @a = split " ";
        if ($chr == $a[1] and $a[2] >= $begin and $a[2] <= $end) {
                if ($a[5] >= 25) {
                        $SNP{$a[2]} = $a[4];
                }
        }
        elsif ($chr < $a[1] or ($chr == $a[1] and $a[2] > $end)) {
                last;
        }
}

close FILE;

print STDERR "got var",(keys %SNP)+0,"\n";

########################################################################

print OUT ">myAcc\n";
for (my $i = $begin; $i <= $end; $i++) {
	
	if (defined($REF{$i})) {
		print OUT $REF{$i};
	}
	elsif (defined($SNP{$i})) {
		print OUT $SNP{$i};
	}
	else {
		print OUT "N";
	}
}
print OUT "\n";

########################################################################


my $database = "ath_eighty";

my $dbh;
&connect_to_db();

my @ecos = ("Agu_1_call","Bak_2_call","Bak_7_call","Cdm_0_call","Del_10_call","Dog_4_call","Don_0_call","Ey15_2_call","Fei_0_call","HKT2_4_call","ICE102_call","ICE104_call","ICE106_call","ICE107_call","ICE111_call","ICE112_call","ICE119_call","ICE120_call","ICE127_call","ICE130_call","ICE134_call","ICE138_call","ICE150_call","ICE152_call","ICE153_call","ICE163_call","ICE169_call","ICE173_call","ICE181_call","ICE1_call","ICE212_call","ICE213_call","ICE216_call","ICE21_call","ICE226_call","ICE228_call","ICE29_call","ICE33_call","ICE36_call","ICE49_call","ICE50_call","ICE60_call","ICE61_call","ICE63_call","ICE70_call","ICE71_call","ICE72_call","ICE73_call","ICE75_call","ICE79_call","ICE7_call","ICE91_call","ICE92_call","ICE93_call","ICE97_call","ICE98_call","Istisu_1_call","Kastel_1_call","Koch_1_call","Lag2_2_call","Leo_1_call","Lerik1_3_call","Memrut_1_call","Mer_6_call","Nie1_2_call","Ped_0_call","Pra_6_call","Qui_0_call","Rue3_1_call","Sha_call","Star_8_call","TueSB30_call","Tuescha9_call","TueV12_call","TueWa1_2_call","Vash_1_call","Vie_0_call","WalhaesB4_call","Xan_1_call","Yeg_1_call");

print "#eco\tcount_ref_matches\tcount_ref_mismatches\tcount_snp_matches\tcount_snp_mismatches\n";

foreach my $eco (@ecos) {

	print STDERR "$eco\n";
	print OUT ">", $eco, "\n";

	my $q = "SELECT chromosome, position, ref, $eco  FROM genome_matrix_calls_25_10_final_idx WHERE chromosome = $chr and position between $begin and $end order by chromosome, position";
	my $sth = $dbh->prepare($q);
	$sth->execute();

	my $count_ref_match = 0;
	my $count_snp_match = 0;
	my $count_ref_mismatch = 0;
	my $count_snp_mismatch = 0;

	while(my $ref = $sth->fetchrow_hashref()) {
	        my $chromosome = $ref->{chromosome};
        	my $position = $ref->{position};
	        my $reference = $ref->{ref};
        	my $allele = $ref->{$eco};

		if ($reference ne $allele) {
			# Case there is a SNP
			if (defined($SNP{$position})) {
				$count_snp_match++;
			}
			elsif (defined($REF{$position})) {
				$count_ref_mismatch++;
			}
		}
		else {
			# Case there is no SNP
			if (defined($REF{$position})) {
                                $count_ref_match++;
                        }
                        elsif (defined($SNP{$position})) {
				$count_snp_mismatch++;
                        }
		}

		print OUT $allele;

	}

	print $eco, "\t", $count_ref_match, "\t", $count_ref_mismatch, "\t", $count_snp_match, "\t", $count_snp_mismatch, "\n"; 
	print OUT "\n";
}







sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) || db_connect_error();
}


