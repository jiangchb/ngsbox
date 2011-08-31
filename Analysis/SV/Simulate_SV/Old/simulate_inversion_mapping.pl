#! /usr/bin/perl

use strict;

my $usage = "simulate_inversion_mapping.pl num";

#my $file = "test_ref_seq.fa"; #shift or die $usage;
#my $index = "test_ref_seq.idx"; 
#my $meta = "test_ref_seq.mta";

my $file = "random_ref_seq.fa"; # 20000 bp
my $index = "random_ref_seq.idx";
my $meta = "random_ref_seq.mta";

my $insert = "pe_insert_freq.txt"; 
my $num = shift or die $usage; 
my $max = shift or die $usage;
my $min = shift or die $usage;
my $offset = 5000; 
my $max_distance_recall = 50; # depends on pe / mp
my $coverage = 20;
my $read_length = 36;

system("rm -rf samples reads mappings sv_pred");

mkdir("samples");
mkdir("reads");
mkdir("mappings");
mkdir("sv_pred");

# Read in reference sequence
open REF, $file or die "file not open";
my $reference = "";

while (my $line = <REF>) {
	chomp($line);
	if (substr($line, 0, 1) ne ">") {
		$reference .= $line;
	}
}
close REF;

# Read in insert dist and create insert probability
open INSERT, $insert or die "insert dist file not open\n";
my %DIST = ();
my $insert_num = 0;
my $insert_min = 1000000000;
my $insert_max = 0;
while(my $line = <INSERT>) {
	my @a = split " ", $line;
	$insert_num += $a[1];
	$DIST{$a[0]} = $a[1];
	if ($a[0] > $insert_max) { $insert_max = $a[0]; }
	if ($a[0] < $insert_min) { $insert_min = $a[0]; }
}
close INSERT;


# Create deletions
open STAT, "> inversions.statistics.".$min."-".$max.".txt";
open LOG, "> inversions.".$min."-".$max.".txt";
open INSERTSIZES, "> inversions.insertsizes.".$min."-".$max.".txt";
open DISCORDANT_INSERTSIZES, "> inversions.insertsizes_discordant.".$min."-".$max.".txt";

my $true_positive = 0;
my $false_negative = 0;
my $false_positive = 0;

for (my $i = 1; $i <= $num; $i++) {

	open SAMPLE, "> samples/sample.".$i.".fa" or die "Cannot open sample file\n";

	######################################################################################################################################################
	# Create sample sequence
	my $start = int(rand(3000)) + $offset; 
	my $length = int(rand($max-$min))+$min;  
	my %nucl_hash=("T" => "A", "G" => "C", "A" => "T", "C" => "G");
	my $inverted_seq = "";
	for (my $nucl = ($start+$length); $nucl >= $start; $nucl--) {
		$inverted_seq .= $nucl_hash{substr($reference, $nucl, 1)};
	}
	my $sample_seq = substr($reference, 0, $start).$inverted_seq.substr($reference, $start+$length, length($reference)-($start+$length));
	print SAMPLE $sample_seq, "\n";
	close SAMPLE;
	
	######################################################################################################################################################
	# Create read pairs
	open READS, "> reads/reads.".$i.".fa";
	open READLOC, "> reads/reads_loc.".$i.".txt";
	my $num_reads = int((length($sample_seq)*$coverage) / $read_length);

	my %READINSERTS = ();
	for (my $i = 0; $i < $num_reads / 2; $i++) {
		my $first_start = (int(rand(length($sample_seq)-1000)));
		my $rand = int(rand($insert_num));
		my $currinsertsize = $insert_min-1;
		while ($rand > 0) {
			$currinsertsize++;
			$rand -= $DIST{$currinsertsize};
		}
		print INSERTSIZES $currinsertsize, "\n";
		$READINSERTS{$i} = $currinsertsize;
		my $second_start = $currinsertsize + $first_start - $read_length;
			
		my $first_sequence = substr($sample_seq, $first_start, $read_length);
		my $second_sequence = reverse(substr($sample_seq, $second_start, $read_length));
		for (my $k=0; $k<length($second_sequence); $k++) {
			if (substr($second_sequence, $k, 1) eq "A") { substr($second_sequence, $k, 1) = "T"; }
			elsif (substr($second_sequence, $k, 1) eq "T") { substr($second_sequence, $k, 1) = "A"; }
			elsif (substr($second_sequence, $k, 1) eq "G") { substr($second_sequence, $k, 1) = "C"; }
			elsif (substr($second_sequence, $k, 1) eq "C") { substr($second_sequence, $k, 1) = "G"; }
		}

		print READS $i, "\t", $first_sequence, "\t", 4, "\t", "qual1", "\t", "qual2", "\t", "qual3", "\n";
		print READS $i, "\t", $second_sequence, "\t", 7, "\t", "qual1", "\t", "qual2", "\t", "qual3", "\n";
		print READLOC $i, "\t", $first_start, "\t", $second_start, "\t", $second_start+$read_length-$first_start+1, "\t", $currinsertsize, "\n";
	}
	close READS;
	close READLOC;	

	# start = 1.deletierte Base, ($start+$length-1) = letzte deletierte Base (Starting from 1!)
	print LOG $i, "\t", $start, "\t", $start+$length-1, "\t", $length;

	#######################################################################################################################################################
	# Map the reads
	system("/ebio/abt6_projects/backup/solexa_tools/genomemapper/genomemapper -i $file -x $index -t $meta -q reads/reads.$i.fa -o mappings/map.$i.list.tmp -E 4 -M 4 -G 3 -e -g 1");
	system("sort -k1n -k2n mappings/map.$i.list.tmp > mappings/map.$i.list");
	system("rm mappings/map.$i.list.tmp");
	#system("perl ~/pgsp/Plot/Mapping/make_plotting_idx.pl --file=mappings/map.$i.list --scale=1000");
	
	# Call structure variants
	system("/ebio/abt6/korbinian/shore/shore structure -i mappings/map.$i.list -o sv_pred/sv.$i.out -s 190 -n 10 -g 5 -b pe_insert_freq.txt -j 20");

	######################################################################################################################################################
	# Parse output


	open OUTPUT, "sv_pred/sv.$i.out/SHORE_SV_inversion.pe.txt";
	my $count_shore_pred = 0; 
	my $pred_start = $max_distance_recall+1;  
	my $pred_end = $max_distance_recall+1; 
	my $found = 0;
	while (my $line = <OUTPUT>) {
		my @a = split " ", $line;
		if (abs($pred_start - $start) > abs($a[1] - $start)) {
			$pred_start = $a[1];
			$pred_end = $a[9];
			$found = 1;
		}
		$count_shore_pred++;
	}
	close OUTPUT;

	if ($found == 1) {
		$true_positive++;
		$false_positive += ($count_shore_pred - 1);
	}
	elsif ($count_shore_pred > 0) {
		$false_positive += ($count_shore_pred - 1);
	}
	elsif ($found == 0) {
		$false_negative++;
	}
	

	# Parse reads in sv cluster 2 get insert dist
	open OUTPUT, "sv_pred/sv.$i.out/readidsincluster";
	while (my $l = <OUTPUT>) {
		my @b = split " ", $l;
		if (defined($READINSERTS{$b[0]})) {
			print DISCORDANT_INSERTSIZES $READINSERTS{$b[0]}, "\n";
		}
	}
	
	# Also SHORE describes the inner positions of the deletions. (Also 1 indexed)
	print LOG "\t", $pred_start, "\t", $pred_end, "\n";
	close OUTPUT;

}
close DISCORDANT_INSERTSIZES;
close INSERTSIZES;
close LOG;

print STAT "Attempts:\t", $num, "\n";
print STAT "TP:\t", $true_positive, "\n";
print STAT "FN:\t", $false_negative, "\n";
print STAT "FP:\t", $false_positive, "\n";

close STAT;

