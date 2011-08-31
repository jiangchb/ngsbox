#! /usr/bin/perl
use warnings;
use strict;

my $usage = "$0 map.list regions readlength \n";
my $map_list = shift or die $usage;
my $regions = shift or die $usage;
my $read_length = shift or die $usage;

if(! defined $read_length) { $read_length = 0; }

### Read target regions
my %REGIONS = ();
open REG, $regions or die "cannot open $regions\n";
while(<REG>) {
	my @a = split " ", $_;
	for (my $i = $a[1]-$read_length; $i<=$a[2]; $i++) {
		$REGIONS{$a[0]."#".$i} = 1;
	}
}
close REG;

### Get reads from target region
open FILE, $map_list or die "cannnot open $map_list\n";
while (<FILE>) {
	my @a = split " ", $_;
	if (defined($REGIONS{$a[0]."#".$a[1]})) {
		print $_;

		### Clean alignment from brackets and read-base
		my $seq = "";
		for (my $i = 0; $i<length($a[2]); $i++) {
			if (substr($a[2], $i, 1) eq "[") {
				$seq .= substr($a[2], $i+2, 1);
				$i+=3;
			}
			else {
				$seq .= substr($a[2], $i, 1);
			}
		}
		print ">$a[3]\n$seq\n";
	}
}


exit(0);
