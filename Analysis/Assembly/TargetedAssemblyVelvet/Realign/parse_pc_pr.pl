#! /usr/bin/perl

my $usage = "perl parse_pc_pr.pl table_file\n";
open FILE, $ARGV[0] or die $usage;

my %CHR_SIZE = ();
$CHR_SIZE{1} = 30432563;
$CHR_SIZE{2} = 19705359;
$CHR_SIZE{3} = 23470805;
$CHR_SIZE{4} = 18585042;
$CHR_SIZE{5} = 26992728;

my $chr;
my $start_real_pcpr;
my $end_real_pcpr;
my $mid_real_pcpr;
my $start_reads;
my $end_reads;
my $read_string_all;
my $read_string1;
my $read_string2;
my %READ_ID = ();

while (my $line = <FILE>) {

	my @a = split " ", $line;

	if ($a[1] < 1) {
		$a[1] = 1;
	}
	if ($a[2] > $CHR_SIZE{$a[0]}) {
		$a[2] = $CHR_SIZE{$a[0]};
	}

	if (defined($chr) and ($chr != $a[0] or $start != $a[1] or $end != $a[2])) {
		print $chr, "\t", $start_reads, "\t", $end_reads + 35, "\t", $start_real_pcpr, "\t", $end_real_pcpr, "\t", $read_string_all, "\t", $read_string1, "\t", $read_string2, "\n";
		$read_string_all = "";
		$read_string1 = "";
		$read_string2 = "";
		$mid_real_pcpr = $a[6] + ($a[7] - $a[6]);
		$start_reads = $a[4];
		$start_real_pcpr = $a[6];
		$end_real_pcpr = $a[7];
	}
	elsif (not defined($start_reads))  {
		$start_reads = $a[4];
		$mid_real_pcpr = $a[6] + ($a[7] - $a[6]);
		$start_real_pcpr = $a[6];
                $end_real_pcpr = $a[7];
	}

	$chr = $a[0];
	$start = $a[1];
	$end = $a[2];
	$end_reads = $a[4];


	$read_string_all .= "," if ($read_string_all ne "");
        $read_string_all .= $a[5];

	my $curr_read_start = $a[4];

	if ($curr_read_start < $mid_real_pcpr - 35) {
		$read_string1 .= "," if ($read_string1 ne "");
		$read_string1 .= $a[5];
	}
	elsif ($curr_read_start > $mid_real_pcpr) {
		$read_string2 .= "," if ($read_string2 ne "");
                $read_string2 .= $a[5];
	}

}

print $chr, "\t", $start_reads, "\t", $end_reads + 35, "\t", $start_real_pcpr, "\t", $end_real_pcpr, "\t", $read_string_all, "\t", $read_string1, "\t", $read_string2, "\n";



