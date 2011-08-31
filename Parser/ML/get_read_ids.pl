#! /usr/bin/perl
use warnings;
use strict;

my $map_list = shift;
my $id_list = shift;

### Store IDs
my %IDS = ();
open IDLIST, $id_list or die "cannot open $id_list\n";
while(<IDLIST>) {
	chomp;
	$IDS{$_} = 1;
}
close IDLIST;

### Get alignments with specified IDs
open FILE, $map_list or die "cannnot open $map_list\n";
my $count = 0;
while (<FILE>) {
	my @a = split " ", $_;

	### DEBUG #######################
	if($count > 100000) {
		print STDERR "$a[0]\t$a[1]\n";
		$count = 0;
	}
	$count++;
	#################################
	
	if ( exists $IDS{$a[3]} ) {
		print $_;
	}
}

exit(0);
