#!/usr/bin/perl

### Convert Shore read files into fastq format
### written by Korbinian Schneeberger, Stephan Ossowski

use strict;
use Getopt::Long;

my %CMD;
my $fl = "";

GetCom();

open FL, $fl or die "Cannot open fl file\n";

while (my $line = <FL>) {
	chomp();
	my @a = split " ", $line;

	print "@".$a[0]."\n";
	print $a[1]."\n";
	print "+\n";
	print $a[4]. "\n";
}



exit(0);

sub GetCom{

  my @usage = ("Usage: $0\n

required:
--fl\tfl formatted file

Will be converted into a fastq file.
	");

	die @usage if ($ARGV[0] eq "");
	GetOptions(\%CMD, "fl=s");

	die("Please specify fl file\n") unless defined($CMD{fl});
  
	$fl = $CMD{fl};
}
	      

