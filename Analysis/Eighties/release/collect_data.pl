#! /usr/bin/perl
use strict;

my $usage = "$0\n";

my $genome_matrix = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2/genome_matrix_new/genome_matrix_max_qual.25.10.wSVdel.wUnseq.txt";
my $sv_base_name = "SV_lib1_version4_final";

if ($genome_matrix eq "") { die "Need to set genome matrix\n"; }

my $base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2";

my $output_base = "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2/genome_matrix_new/release_2009_12_27";
if (not (-e $output_base)) {
	mkdir($output_base);
}

#my @ECOTYPES = ("Agu-1","Bak-2");
#my @ECOTYPES = ("Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE173","ICE181","ICE1","ICE212","ICE213","ICE216","ICE21","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE7","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");
my @ECOTYPES = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE173","ICE181","ICE1","ICE212","ICE213","ICE216","ICE21","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE7","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");

print STDERR "In total:", @ECOTYPES+0, "\n";

my $cp_insertions = 1;
my $cp_unsequenced = 1;
my $cp_quality_calls = 1;
my $cp_sv_calls = 1;
my $parse_genome_matrix = 1;


for (my $i = 0; $i <@ECOTYPES; $i++) {
	print STDERR "Analysing ".$ECOTYPES[$i]."\n";

	####################################################
	# Create ecotype folder
	####################################################
	my $output_ecotype = $output_base."/".$ECOTYPES[$i];
	if (not (-e $output_ecotype)) {
		mkdir $output_ecotype;
	}

	####################################################
	# Copy short insertions
	####################################################
	my $analysis = $base."/".$ECOTYPES[$i]."/AlignmentFolder/Analysis_q01/ConsensusAnalysis/";
	my $short_ins = $analysis."/insertion.txt";
	my $ins_cmd = "cp -l $short_ins $output_ecotype/.";
	if ($cp_insertions == 1) {
		system($ins_cmd);
	}

	####################################################
	# Copy unsequenced
	####################################################
	my $unseq = $analysis."/unsequenced.txt";
	my $unseq_cmd = "cp -l $unseq $output_ecotype/.";
	if ($cp_unsequenced == 1) {
		system($unseq_cmd);
	}

	####################################################
	# Copy quality_variants.txt and quality_reference.txt
	# Read in files
	####################################################
	my $ref = $analysis."/quality_reference.txt";
        my $ref_cmd = "cp -l $ref $output_ecotype/.";
	if ($cp_quality_calls == 1) {
       		system($ref_cmd);		
	}
	my $var = $analysis."/quality_variant.txt";
        my $var_cmd = "cp -l $var $output_ecotype/.";
	if ($cp_quality_calls == 1) {
        	system($var_cmd);
	}

	my %REF = ();
	if ($parse_genome_matrix == 1) {
		open REFFILE, $ref or die "Cannot open reference file: ".$ref."\n";
		while (my $line = <REFFILE>) {
			my @a = split " ", $line;
			$REF{$a[1]."#".$a[2]} = $line;	
		}
		close REFFILE;
	}
	
	my %VAR = ();
	if ($parse_genome_matrix == 1) {
        	open VARFILE, $var or die "Cannot open variant file: ".$var."\n";
 	       	while (my $line = <VARFILE>) {
        	        my @a = split " ", $line;
                	$VAR{$a[1]."#".$a[2]} = $line;
		}
        	close VARFILE;
        }


	####################################################
	# Copy SV
	####################################################
	my $sv = $base."/".$ECOTYPES[$i]."/AlignmentFolder/$sv_base_name/";
	my $sv_del = $sv."/SV_deletion_high_quality.lib1.PE.txt";
	my $sv_ins = $sv."/SV_insertion_high_quality.lib1.PE.txt";
	my $sv_inv = $sv."/SV_inversion_high_quality.lib1.PE.txt";
	my $sv_del_cmd = "cp -l $sv_del $output_ecotype/.";
	my $sv_ins_cmd = "cp -l $sv_ins $output_ecotype/.";
	my $sv_inv_cmd = "cp -l $sv_inv $output_ecotype/.";
	if (-e $sv_del) {
		system($sv_del_cmd);
	}
	if (-e $sv_ins) {
		system($sv_ins_cmd);
	}
	if (-e $sv_inv) {
		system($sv_inv_cmd);
	}

	####################################################
	# Parse genome matrix
	####################################################
	if ($parse_genome_matrix == 1) {
		open GMFILE, $genome_matrix or die "Cannot open genome matrix\n";
		open SNPFILE, ">".$output_ecotype."/quality_variant_genome_matrix.txt";
		open REFFILE, ">".$output_ecotype."/quality_reference_genome_matrix.txt";

		my @Z = ();
		my @ZN = ();

		while (my $line = <GMFILE>) {
			my @a = split " ", $line;
			if ($a[2] eq $a[$i+3]) {
				print REFFILE $REF{$a[0]."#".$a[1]};
			}
			else {
				if ($a[$i+3] eq "A" || $a[$i+3] eq "C" || $a[$i+3] eq "G" || $a[$i+3] eq "T" || $a[$i+3] eq "-") {
					print SNPFILE $VAR{$a[0]."#".$a[1]};
				}
				elsif ($a[$i+3] eq "Z" || $a[$i+3] eq "N") {
					push @ZN, $a[0]."#".$a[1];
					if ($a[$i+3] eq "Z") {
						push @Z, $a[0]."#".$a[1];
					}
				}
			}
		}
	
		close GMFILE;	
		close SNPFILE;	
		close REFFILE;	
	

		####################################################
		# Write highly_diverged_regions.txt and inaccessible_regions.txt
		####################################################
		my $hdr = ">".$output_ecotype."/highly_diverged_regions.txt";		
		my $noa = ">".$output_ecotype."/inaccessible_regions.txt";	

		parse_array(\@Z, $hdr);
		parse_array(\@ZN, $noa);
	}
	
}

sub parse_array {
	my ($array, $file_name) = @_;
	open OUT, ">".$file_name;

	my $last_chr = -1;
	my $last_pos = -1;
	my $start = 0;

	for (my $i = 0; $i < @$array; $i++) {
		my ($chr, $pos) = split "#", ${@$array}[$i];
		if ($last_chr != $chr) {
			if ($start != 0) {
				print OUT $last_chr, "\t", $start, "\t", $last_pos, "\n";
			}
			$start = $pos;
		}	
		elsif ($last_pos+1 != $pos) {
			if ($start != 0) {
				print OUT $last_chr, "\t", $start, "\t", $last_pos, "\n";
			}
			$start = $pos;
		}

		$last_chr = $chr;
		$last_pos = $pos;
	}

	close OUT;	
}















