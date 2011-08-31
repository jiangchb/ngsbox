#! /usr/bin/perl
use strict;

my @econames = ("Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE173","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Memrut-1","Mer-6","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1","Sha","Star-8","TueSB30","Tuescha9","TueV12","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");

my %NAME2ID = ();
my %ID2NAME = ();

$NAME2ID{"Agu-1"} = 0; $ID2NAME{0} = "Agu-1"; $NAME2ID{"Bak-2"} = 1; $ID2NAME{1} = "Bak-2"; $NAME2ID{"Bak-7"} = 2; $ID2NAME{2} = "Bak-7"; $NAME2ID{"Cdm-0"} = 3; $ID2NAME{3} = "Cdm-0"; $NAME2ID{"Del-10"} = 4; $ID2NAME{4} = "Del-10"; $NAME2ID{"Dog-4"} = 5; $ID2NAME{5} = "Dog-4"; $NAME2ID{"Don-0"} = 6; $ID2NAME{6} = "Don-0"; $NAME2ID{"Ey15-2"} = 7; $ID2NAME{7} = "Ey15-2"; $NAME2ID{"Fei-0"} = 8; $ID2NAME{8} = "Fei-0"; $NAME2ID{"HKT2.4"} = 9; $ID2NAME{9} = "HKT2.4"; $NAME2ID{"ICE1"} = 10; $ID2NAME{10} = "ICE1"; $NAME2ID{"ICE102"} = 11; $ID2NAME{11} = "ICE102"; $NAME2ID{"ICE104"} = 12; $ID2NAME{12} = "ICE104"; $NAME2ID{"ICE106"} = 13; $ID2NAME{13} = "ICE106"; $NAME2ID{"ICE107"} = 14; $ID2NAME{14} = "ICE107"; $NAME2ID{"ICE111"} = 15; $ID2NAME{15} = "ICE111"; $NAME2ID{"ICE112"} = 16; $ID2NAME{16} = "ICE112"; $NAME2ID{"ICE119"} = 17; $ID2NAME{17} = "ICE119"; $NAME2ID{"ICE120"} = 18; $ID2NAME{18} = "ICE120"; $NAME2ID{"ICE127"} = 19; $ID2NAME{19} = "ICE127"; $NAME2ID{"ICE130"} = 20; $ID2NAME{20} = "ICE130"; $NAME2ID{"ICE173"} = 21; $ID2NAME{21} = "ICE173"; $NAME2ID{"ICE134"} = 22; $ID2NAME{22} = "ICE134"; $NAME2ID{"ICE138"} = 23; $ID2NAME{23} = "ICE138"; $NAME2ID{"ICE150"} = 24; $ID2NAME{24} = "ICE150"; $NAME2ID{"ICE152"} = 25; $ID2NAME{25} = "ICE152"; $NAME2ID{"ICE153"} = 26; $ID2NAME{26} = "ICE153"; $NAME2ID{"ICE163"} = 27; $ID2NAME{27} = "ICE163"; $NAME2ID{"ICE169"} = 28; $ID2NAME{28} = "ICE169"; $NAME2ID{"ICE181"} = 29; $ID2NAME{29} = "ICE181"; $NAME2ID{"ICE21"} = 30; $ID2NAME{30} = "ICE21"; $NAME2ID{"ICE212"} = 31; $ID2NAME{31} = "ICE212"; $NAME2ID{"ICE213"} = 32; $ID2NAME{32} = "ICE213"; $NAME2ID{"ICE216"} = 33; $ID2NAME{33} = "ICE216"; $NAME2ID{"ICE226"} = 34; $ID2NAME{34} = "ICE226"; $NAME2ID{"ICE228"} = 35; $ID2NAME{35} = "ICE228"; $NAME2ID{"ICE29"} = 36; $ID2NAME{36} = "ICE29"; $NAME2ID{"ICE33"} = 37; $ID2NAME{37} = "ICE33"; $NAME2ID{"ICE36"} = 38; $ID2NAME{38} = "ICE36"; $NAME2ID{"ICE49"} = 39; $ID2NAME{39} = "ICE49"; $NAME2ID{"ICE50"} = 40; $ID2NAME{40} = "ICE50"; $NAME2ID{"ICE60"} = 41; $ID2NAME{41} = "ICE60"; $NAME2ID{"ICE61"} = 42; $ID2NAME{42} = "ICE61"; $NAME2ID{"ICE63"} = 43; $ID2NAME{43} = "ICE63"; $NAME2ID{"ICE7"} = 44; $ID2NAME{44} = "ICE7"; $NAME2ID{"ICE70"} = 45; $ID2NAME{45} = "ICE70"; $NAME2ID{"ICE71"} = 46; $ID2NAME{46} = "ICE71"; $NAME2ID{"ICE72"} = 47; $ID2NAME{47} = "ICE72"; $NAME2ID{"ICE73"} = 48; $ID2NAME{48} = "ICE73"; $NAME2ID{"ICE75"} = 49; $ID2NAME{49} = "ICE75"; $NAME2ID{"ICE79"} = 50; $ID2NAME{50} = "ICE79"; $NAME2ID{"ICE91"} = 51; $ID2NAME{51} = "ICE91"; $NAME2ID{"ICE92"} = 52; $ID2NAME{52} = "ICE92"; $NAME2ID{"ICE93"} = 53; $ID2NAME{53} = "ICE93"; $NAME2ID{"ICE97"} = 54; $ID2NAME{54} = "ICE97"; $NAME2ID{"ICE98"} = 55; $ID2NAME{55} = "ICE98"; $NAME2ID{"Istisu-1"} = 56; $ID2NAME{56} = "Istisu-1"; $NAME2ID{"Kastel-1"} = 57; $ID2NAME{57} = "Kastel-1"; $NAME2ID{"Koch-1"} = 58; $ID2NAME{58} = "Koch-1"; $NAME2ID{"Lag2.2"} = 59; $ID2NAME{59} = "Lag2.2"; $NAME2ID{"Leo-1"} = 60; $ID2NAME{60} = "Leo-1"; $NAME2ID{"Lerik1-3"} = 61; $ID2NAME{61} = "Lerik1-3"; $NAME2ID{"Memrut-1"} = 62; $ID2NAME{62} = "Memrut-1"; $NAME2ID{"Mer-6"} = 63; $ID2NAME{63} = "Mer-6"; $NAME2ID{"Nie1-2"} = 64; $ID2NAME{64} = "Nie1-2"; $NAME2ID{"Ped-0"} = 65; $ID2NAME{65} = "Ped-0"; $NAME2ID{"Pra-6"} = 66;$ID2NAME{66} = "Pra-6"; $NAME2ID{"Qui-0"} = 67; $ID2NAME{67} = "Qui-0"; $NAME2ID{"Rue3-1"} = 68; $ID2NAME{68} = "Rue3-1"; $NAME2ID{"Sha"} = 69; $ID2NAME{69} = "Sha"; $NAME2ID{"Star-8"} = 70; $ID2NAME{70} = "Star-8"; $NAME2ID{"TueSB30"} = 71; $ID2NAME{71} = "TueSB30"; $NAME2ID{"Tuescha9"} = 72; $ID2NAME{72} = "Tuescha9"; $NAME2ID{"TueV12"} = 73; $ID2NAME{73} = "TueV12"; $NAME2ID{"TueWa1-2"} = 74; $ID2NAME{74} = "TueWa1-2"; $NAME2ID{"Vash-1"} = 75; $ID2NAME{75} = "Vash-1"; $NAME2ID{"Vie-0"} = 76; $ID2NAME{76} = "Vie-0"; $NAME2ID{"WalhaesB4"} = 77; $ID2NAME{77} = "WalhaesB4"; $NAME2ID{"Xan-1"} = 78; $ID2NAME{78} = "Xan-1"; $NAME2ID{"Yeg-1"} = 79; $ID2NAME{79} = "Yeg-1";

### Set up global

my $ecostring_raw = "";
for (my $i = 0; $i < 80; $i++) {
	$ecostring_raw .= "0";
}

my $sv_extension_blocking = 500;
my $sv_size_blocking = 50;

my $usage = "$0 SVfrequency\n";
my $file = shift or die $usage;

### Read in SVs by Frequency and Pattern

my %ID_SV = ();
my %SV = ();
my %BLOCK = ();

open FILE, $file or die "Cannot open file\n";
while (my $line = <FILE>) {
	my @a = split " ", $line;
	$ID_SV{$a[11]} .= $line;
}
close FILE;

print STDERR "reads sv file\n";

foreach my $id (keys %ID_SV) {
	my @a = split "\n", $ID_SV{$id};
	my $chr;
	my $start=-1;
	my $end=-1;
	my $ecostring = $ecostring_raw; 
	for (my $i = 0; $i < @a; $i++) {
		# store accession
		my @b = split " ", $a[$i];
		substr($ecostring, $NAME2ID{$b[0]}, 1, "1"); 
		# store location
		$chr = $b[3];
		if ($start == -1 || $start > $b[4]) {
			$start = $b[4];
		}
		if ($end == -1 || $end < $b[5]) {
			$end = $b[5];
		}
		# Set up blocked regions
		if ($end-$start+1 >= $sv_size_blocking) {
			for (my $j = $b[4]-$sv_extension_blocking; $j <= $b[5]+$sv_extension_blocking; $j++) {
				$BLOCK{$b[3]."#".$j} = 1;
			}
		}
	}
	#my $ecostring = join ",", sort @ecos;
	my $chrpos = ($chr*100000000)+$start;
	my $chrend = ($chr*100000000)+$end;
	$SV{@a+0}{$ecostring} .= $chrpos."-".$chrend."#";
}

print STDERR "clac sv freq.\n";
print STDERR "blocked regions: ", (keys %BLOCK)+0, "\n";

#open OUT, ">tmp.out";
#foreach my $key (keys %SV) {
#	foreach my $key2 (keys %{$SV{$key}}) {
#		print OUT $key, "\t", $key2, "\t", $SV{$key}{$key2}, "\n";
#	}
#}
#exit(1);

### Parse genome matrix for identical patterns

open FILE, "/ebio/abt6_projects2/backup/data/solexa_analysis/ATH/Genome/EIGHTY2/genome_matrix_new/genome_matrix_max_qual.25.10.wSVdel.wUnseq.txt";
#open FILE, "test_input_gm";
my %SNP = ();
my $count = 0;
while (my $line = <FILE>) {
	$count++;

	my ($chr, $pos) = split " ", $line;
	print STDERR $chr, "\t", $pos, "\n" if $count%100000==0;

	if (not defined($BLOCK{$chr."#".$pos})) {
		my @a = split " ", $line;
		my $ecostring = $ecostring_raw;
		my $count = 0;
		for (my $i = 3; $i < @a; $i++) {
			if ($a[2] ne $a[$i] and ($a[$i] eq "A" || $a[$i] eq "C" || $a[$i] eq "G" || $a[$i] eq "T")) {
				substr($ecostring, $i-3, 1, "1");
				$count++;
			}
		}
		my $chrpos = ($a[0]*100000000)+$a[1];
		$SNP{$count}{$ecostring} .= $chrpos."#"; 
	}
}

print STDERR "read genome matrix\n";

### 
open OUT, ">SV_SNP_pattern.e".$sv_extension_blocking.".s".$sv_size_blocking.".txt";
#open OUT, ">SV_SNP_pattern.e".$sv_extension_blocking.".s".$sv_size_blocking.".test";

for (my $freq = 1; $freq <= 80; $freq++) {
	foreach my $ecostring (keys %{$SV{$freq}}) {
		my $match_out = "";
		my $fuzzy_out = "";

		if (defined($SNP{$freq}{$ecostring})) {
			$match_out = $SNP{$freq}{$ecostring}, "\n";
		}

		$fuzzy_out = get_fuzzy_matches($ecostring, $freq);
		
		# Print if results:
		if ($match_out ne "" or $fuzzy_out ne "") {
			print OUT $freq, "\t", get_econames($ecostring, $freq), "\t", $ecostring, "\n";
			print OUT $SV{$freq}{$ecostring}, "\n";
			print OUT $match_out, "\n";
			print OUT $fuzzy_out, "\n";	
		}

	}
}

sub get_fuzzy_matches {
        my ($bin, $freq) = @_;

	my $str = "";

        for (my $i = 0; $i < length($bin); $i++) {
		my $tmp_bin = $bin;
		my $tmp_freq = $freq;
		if (substr($tmp_bin, $i, 1) eq "0") {
			substr($tmp_bin, $i, 1, "1");
			$tmp_freq++;
		}
		else { 
			substr($tmp_bin, $i, 1, "0");
			$tmp_freq--;
		}
		if (defined($SNP{$tmp_freq}{$tmp_bin})) {
			if (int(rand(100)) == 0) {
				$str .= $SNP{$tmp_freq}{$tmp_bin}."#";
			}
		}
	}

	return $str;
}

sub get_econames {
	my ($bin, $freq) = @_;
	
	my $str = "";

	for (my $i = 0; $i < length($bin); $i++) {
		if (substr($bin, $i, 1) eq "1") {
			$str .= "," if $str ne "";
			$str .= $ID2NAME{$i};
		}
	}

	return $str;
}






