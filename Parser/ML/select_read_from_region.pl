#! /usr/bin/perl
use warnings;
use strict;

my $map_list = shift;
my $chr      = shift;
my $start    = shift;
my $end      = shift;
my $type     = shift; # fasta oder fastq
 

### Get reads from target region
open FILE, $map_list or die "cannnot open $map_list\n";
my $written = 0;
my %READS = ();
while (<FILE>) {
	my @a = split " ", $_;

	if( $a[0] == $chr && $a[1] >= $start && $a[1] <= $end) {
		if (not defined($READS{$a[3].".".$a[9]})) {
			if ($type eq "fasta") {
				print ">", $a[3], "\n", get_seq($a[2]), "\n";
			}
			elsif ($type eq "fastq") {
				print "@", $a[3], "\n", get_seq($a[2]), "\n+\n", $a[11], "\n";
			}
			$READS{$a[3].".".$a[9]} = 1;
		}
		$written = 1;
	}
	else {
		if ($written == 1) {
			exit(0);
		}
	}
}


exit(0);

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

