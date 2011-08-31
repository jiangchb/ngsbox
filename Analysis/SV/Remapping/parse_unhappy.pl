#! /usr/bin/perl

use strict;

my $usage = "\n$0 maplist lib folder\n\n";
my $file = shift or die $usage;
my $lib = shift or die $usage;
my $base = shift or die $usage;

mkdir "$base";
mkdir "$base/1";
mkdir "$base/1/1";
mkdir "$base/1/2";

my %LEN_USED_1 = ();
my %LEN_USED_2 = ();
my %LEN_BATCH_1 = ();
my %LEN_BATCH_2 = ();

my %READS = ();

open FILE, $file or die $usage;

my $curr_chr = -1;

my $num1 = 0;
my $num2 = 0;

while (<FILE>) {
	my @a = split " ";
	if ($curr_chr != $a[0]) {
		%READS = ();
		$curr_chr = $a[0];
	}

	if ($a[1] % 1000000 == 0) {
		print STDERR $a[0], "\t", $a[1], "\n";
	}

	if ($a[6] == 1 and (($lib == 1 and $a[9] == 4 || $a[9] == 7) or ($lib == 2 and ($a[9] == 10 || $a[9] == 13)) or ($lib == 3 and ($a[9] == 16 || $a[9] == 19)))) {
		if (defined($READS{$a[3]})) {

			my $seq1 = $READS{$a[3]};
			my $seq2 = "";
			if ($a[4] eq "D") { 
				$seq2 = get_seq($a[2]);
			}
			else {
				$seq2 = rev_comp(get_seq($a[2]));
			}

			$num1++; 
			$num2++;
			
			#my $len1 = int($num1 / 100000) + 1;
			#my $len2 = int($num2 / 100000) + 1;
			
			my $len1 = length($seq1);
			my $len2 = length($seq2);

			if (($a[9] >= 3 and $a[9] <= 5) or ($a[9] >= 9 and $a[9] <= 11)) {
	
				if (not defined($LEN_USED_1{$len2})) {
					$LEN_USED_1{$len2} = 1;
					$LEN_BATCH_1{$len2} = 0; 
				}
				else {
					$LEN_USED_1{$len2}++;
					if ($LEN_USED_1{$len2}%50000 == 0) {
						$LEN_USED_1{$len2} = 1;
						$LEN_BATCH_1{$len2}++;
					}
				}
				my $batch1 = $LEN_BATCH_1{$len2};
				my $folder1 = "$base/1/1/length_".$len2."_$batch1";
				if (not -e $folder1) {
					mkdir $folder1;
				}
	
				if (not defined($LEN_USED_2{$len1})) {
					$LEN_USED_2{$len1} = 1;
					$LEN_BATCH_2{$len1} = 0;
                                }
				else {
					$LEN_USED_2{$len1}++;
                                        if ($LEN_USED_2{$len1}%50000 == 0) {
                                                $LEN_USED_2{$len1} = 1;
                                                $LEN_BATCH_2{$len1}++;
                                        }
				}
				my $batch2 = $LEN_BATCH_2{$len1};
                                my $folder2 = "$base/1/2/length_".$len1."_$batch2";
                                if (not -e $folder2) {
                                        mkdir $folder2;
                                }

				open OUT1, ">> $folder1/reads_0.fl";
				open OUT2, ">> $folder2/reads_0.fl";
				print OUT1 $a[3], "\t", $seq2, "\t1\tq1\tq2\tq3\n";
				print OUT2 $a[3], "\t", $seq1, "\t2\tq1\tq2\tq3\n";
				close OUT1;
				close OUT2;
			}
			else {
				if (not defined($LEN_USED_1{$len1})) {
                                        $LEN_USED_1{$len1} = 1;
                                        $LEN_BATCH_1{$len1} = 0;
                                }
                                else {
                                        $LEN_USED_1{$len1}++;
                                        if ($LEN_USED_1{$len1}%50000 == 0) {
                                                $LEN_USED_1{$len1} = 1;
                                                $LEN_BATCH_1{$len1}++;
                                        }
                                }
                                my $batch1 = $LEN_BATCH_1{$len1};
                                my $folder1 = "$base/1/1/length_".$len1."_$batch1";
                                if (not -e $folder1) {
                                        mkdir $folder1;
                                }

				if (not defined($LEN_USED_2{$len2})) {
                                        $LEN_USED_2{$len2} = 1;
                                        $LEN_BATCH_2{$len2} = 0;
                                }
                                else {
                                        $LEN_USED_2{$len2}++;
                                        if ($LEN_USED_2{$len2}%50000 == 0) {
                                                $LEN_USED_2{$len2} = 1;
                                                $LEN_BATCH_2{$len2}++;
                                        }
                                }
                                my $batch2 = $LEN_BATCH_2{$len2};
                                my $folder2 = "$base/1/2/length_".$len2."_$batch2";
                                if (not -e $folder2) {
                                        mkdir $folder2;
                                }


                                open OUT1, ">> $folder1/reads_0.fl";
                                open OUT2, ">> $folder2/reads_0.fl";
				print OUT1 $a[3], "\t", $seq1, "\t1\tq1\tq2\tq3\n";
                       	        print OUT2 $a[3], "\t", $seq2, "\t2\tq1\tq2\tq3\n";
				close OUT1;
                                close OUT2;
			}

			delete $READS{$a[3]};
		}
		else  {
			if ($a[4] eq "D") {
				$READS{$a[3]} = get_seq($a[2]);
			}
			else {
				$READS{$a[3]} = rev_comp(get_seq($a[2]));
			}
		}	
	}
}

sub get_seq {
	my ($align) = @_;

	my $seq = "";
        for (my $i = 0; $i<length($align); $i++) {
        	if (substr($align, $i, 1) eq "[") {
                	if (substr($align, $i+2, 1) ne "-") {
                        	$seq .= substr($align, $i+2, 1);
                        }
                        $i+=3;
                }
                else {
                	$seq .= substr($align, $i, 1);
                }
        }

	return $seq;

}

sub rev_comp {
        my $seq = shift;
        my $rev = "";

        for (my $i = 0; $i<length($seq); $i++) {
                $rev .= rev_comp_base(substr($seq, $i, 1));
        }

        my $s = reverse $rev;

        return $s;
}

sub rev_comp_base {
        my ($base) = @_;

        $base =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

        return $base;
}



