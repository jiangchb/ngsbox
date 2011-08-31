#! /usr/bin/perl
use strict;
use warnings;

my $usage = "$0 fastafile\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

my %entries = ();
my $id = "";

while (my $line = <FILE>) {

	if (substr($line, 0, 1) eq ">") {
		my @e = split(" ", $line);
		$id = substr($e[0], 1);

		# Human specific modifications
		if($id eq "X")     { $id = 23; }
		elsif($id eq "Y")  { $id = 24; }
		elsif($id eq "MT") { $id = 25; }

		$entries{$id} = $line;
	}
	else {
		$entries{$id} .= $line;
	}
}

foreach my $header ( sort {$a<=>$b} keys %entries ) {
	print $entries{$header};
}
	

