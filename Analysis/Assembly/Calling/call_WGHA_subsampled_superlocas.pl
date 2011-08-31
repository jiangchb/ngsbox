#! /usr/bin/perl
use strict;

# This script assumes that the AssemblyFolder is on the same level as the AlignmentFolder

my $usage = "$0 AssemFolder shorebinary chrset left_arm right_arm\n";
my $assemfolder = shift or die $usage;
my $shore = shift;
my $chrset = shift;
my $left_arm_flag = shift;
my $right_arm_flag = shift;

print "Assemfolder:", $assemfolder, "\n";
print "shorebin:", $shore, "\n";

chdir($assemfolder);


# TAIR8:
my @left_end = (13700000, 2450000, 11300000, 1800000, 11000000);
my @right_start = (15900000, 5500000, 14300000, 5150000, 13350000);
my @right_end = (30432563, 19705359, 23470805, 18585042, 26992728);

my @chrs = split ",", $chrset;

foreach my $chr (@chrs) {

	chdir("chr".$chr);	

	## Left Arm

	if ($left_arm_flag == 1) {

		if (not -e "Assembly_Left_Arm_superlocas_0".$chr) {
			system("~/shore/startshore.sh WGHA -L ../libraries.txt -i ../map.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Left_Arm_superlocas_0$chr -s -e ../left_over/subsam_left_over_1_2.lib1.fq,../left_over/subsam_left_over_1_2.lib2.fq,../left_over/subsam_left_over_1_2.lib3.fq,../left_over/subsam_left_over_single.all.fq -I 1 -U 50 -T 25 -m 12000 -j 21 -M 60 -g 0 -C $chr -S 1 -E ".$left_end[$chr-1]);
		}
	
	}

	## Right Arm

	if ($right_arm_flag == 1) {

		if (not -e "Assembly_Right_Arm_superlocas_0".$chr) {
			system("~/shore/startshore.sh WGHA -L ../libraries.txt -i ../map.list -l ../left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o Assembly_Right_Arm_superlocas_0$chr -s -e ../left_over/subsam_left_over_1_2.lib1.fq,../left_over/subsam_left_over_1_2.lib2.fq,../left_over/subsam_left_over_1_2.lib3.fq,../left_over/subsam_left_over_single.all.fq -I 1 -U 50 -T 25 -m 12000 -j 21 -M 60 -g 0 -C $chr -S ".$right_start[$chr-1]." -E ".$right_end[$chr-1]);
	        }

	}

	chdir("..");

}


