#! /usr/bin/perl
use strict;

# This script assumes that the AssemblyFolder is on the same level as the AlignmentFolder

my $usage = "$0 AssemFolder  velvet abyss  euler  superlocas  shorebinary  kmerset  chrset  left_arm right_arm\n";
my $assemfolder = shift or die $usage;
my $run_velvet = shift;
my $run_abyss = shift;
my $run_euler = shift;
my $run_superlocas = shift;
my $shore = shift;
my $kmerset = shift;
my $chrset = shift;
my $left_arm_flag = shift;
my $right_arm_flag = shift;

print "Assemfolder:", $assemfolder, "\n";
print "velvet:", $run_velvet, "\n";
print "abyss:", $run_abyss, "\n";
print "euler:", $run_euler, "\n";
print "superlocas:", $run_superlocas, "\n";
print "shorebin:", $shore, "\n";
print "kmerset:", $kmerset, "\n";

my $oversampled="../Analysis_q03/ConsensusAnalysis/supplementary_data/oversampled.txt";

chdir($assemfolder);

my $kmers_velvet_abyss = "";
my $kmers_euler = "";

if ($kmerset==1) {
	# Descent read length (Ler-1):
	$kmers_velvet_abyss = "23,27,31,35,39,43,47,51";
	$kmers_euler = "21,22,23,24,25,26,27,28";
}
elsif ($kmerset==2) {
	# Mixed read length: (Kro-0)
	$kmers_velvet_abyss = "19,21,23,27,31,37,41,45";
	$kmers_euler = "17,19,21,23,25,26,27,28";
}
elsif ($kmerset==3) {
	# Short read length: (Bur-0, C24, Col-0)
	$kmers_velvet_abyss = "17,19,21,23,25,27,29,31";
	$kmers_euler = "17,18,20,21,23,24,26,28";
}

elsif ($kmerset==4) {
	$kmers_velvet_abyss = "23,27,31,35";
	$kmers_euler = "17,18,20,21";
}
elsif ($kmerset==5) {
	$kmers_velvet_abyss = "39,43,47,51";
	$kmers_euler = "23,24,26,28";
}
elsif ($kmerset==6) {
        $kmers_velvet_abyss = "23,27";
        $kmers_euler = "17,18";
}
elsif ($kmerset==7) {
	$kmers_velvet_abyss = "23";
	$kmers_euler = "17";
}
else {
	die("Need kmerset set to 1, 2, 3, 4 or 5.\n");
}

# TAIR8:
my @left_end = (13700000, 2450000, 11300000, 1800000, 11000000);
my @right_start = (15900000, 5500000, 14300000, 5150000, 13350000);
my @right_end = (30432563, 19705359, 23470805, 18585042, 26992728);

my @chrs = split ",", $chrset;

foreach my $chr (@chrs) {

	chdir("chr".$chr);	



	### Left Arm
	if ($left_arm_flag == 1) {

		# Velvet
		if (not -e "Assembly_Left_Arm_velvet_0".$chr and $run_velvet == 1) {	
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Left_Arm_velvet_0$chr -v -m 12000 -j $kmers_velvet_abyss -M 52 -g 0 -k 8 -z -C $chr -S 1 -E ".$left_end[($chr-1)]." -X $oversampled");
		}

		# Abyss
		if (not -e "Assembly_Left_Arm_abyss_0".$chr and $run_abyss == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Left_Arm_abyss_0$chr -y -m 12000 -j $kmers_velvet_abyss -M 52 -g 0 -k 8 -z -C $chr -S 1 -E ".$left_end[($chr-1)]." -X $oversampled");
		}

		# Euler-SR
		if (not -e "Assembly_Left_Arm_euler_0".$chr and $run_euler == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Left_Arm_euler_0$chr -u -m 12000 -j $kmers_euler -M 52 -g 0 -k 8 -z -C $chr -S 1 -E ".$left_end[$chr-1]." -X $oversampled");
		}

		# Superlocas
		if (not -e "Assembly_Left_Arm_superlocas_0".$chr and $run_superlocas == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Left_Arm_superlocas_0$chr -s -e ../left_over/left_over_1_2.lib1.fq,../left_over/left_over_1_2.lib2.fq,../left_over/left_over_1_2.lib3.fq -m 12000 -j 21 -M 52 -g 0 -k 8 -z -C $chr -S 1 -E ".$left_end[$chr-1]." -X $oversampled");
		}
	}



	### Right Arm
	if ($right_arm_flag == 1) {

		# Velvet
		if (not -e "Assembly_Right_Arm_velvet_0".$chr and $run_velvet == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Right_Arm_velvet_0$chr -v -m 12000 -j $kmers_velvet_abyss -M 52 -g 0 -k 8 -z -C $chr -S ".$right_start[$chr-1]." -E ".$right_end[$chr-1]." -X $oversampled");
		}

		# Abyss
		if (not -e "Assembly_Right_Arm_abyss_0".$chr and $run_abyss == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Right_Arm_abyss_0$chr -y -m 12000 -j $kmers_velvet_abyss -M 52 -g 0 -k 8 -z -C $chr -S ".$right_start[$chr-1]." -E ".$right_end[$chr-1]." -X $oversampled");
		}

		# Euler-SR
		if (not -e "Assembly_Right_Arm_euler_0".$chr and $run_euler == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Right_Arm_euler_0$chr -u -m 12000 -j $kmers_euler -M 52 -g 0 -k 8 -z -C $chr -S ".$right_start[$chr-1]." -E ".$right_end[$chr-1]." -X $oversampled");
		}

		# Superlocas
		if (not -e "Assembly_Right_Arm_superlocas_0".$chr and $run_superlocas == 1) {
			system("$shore WGHA -L ../libraries.txt -i map_chr$chr.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Right_Arm_superlocas_0$chr -s -e ../left_over/left_over_1_2.lib1.fq,../left_over/left_over_1_2.lib2.fq,../left_over/left_over_1_2.lib3.fq -m 12000 -j 21 -M 52 -g 0 -k 8 -z -C $chr -S ".$right_start[$chr-1]." -E ".$right_end[$chr-1]." -X $oversampled");
	        }
	}


	chdir("..");
}
