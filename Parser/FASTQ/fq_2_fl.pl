#!/usr/bin/perl

# Convert FASTQ to Shore sequence read format
# Written by Korbinian Schneeberger, Stephan Ossowski


use strict;
use Getopt::Long;

my %CMD;
my $fastq = "";
my $length = "";

GetCom();


### prepare fake qcal and chas
my $qcal = ""; 
my $chas = "";
for(my $i = 0; $i < $length; $i++) {
	$qcal .= "Z";
	$chas .= "Z";
}

open FASTQ, $fastq or die "Cannot open fastq file\n";

while (<FASTQ>) {
	chomp;
	my $line = $_;

	if($line =~ /\@/) {
		my @members = split(/_/, $line);

		for (my $i = 3; $i < 6; $i++) {
			$members[$i] =~ s/-//g;
			if( length($members[$i]) == 1 ) { $members[$i] = "000" . $members[$i]; }
			elsif( length($members[$i]) == 2 ) { $members[$i] = "00" . $members[$i]; }
			elsif( length($members[$i]) == 3 ) { $members[$i] = "0" . $members[$i]; }
		}
		my $id = $members[2] . $members[3] . $members[4] . $members[5];

		my $seq = <FASTQ>;
		chomp($seq);
		print "$id\t$seq\t";
	}
	elsif($line =~ /^+/) {
		my $qval = <FASTQ>;
		chomp($qval);
		print "$qval\t$qcal\t$chas\n";
	}
	else { print "file format not correct!\n"; exit(0); }
}




exit(0);

sub GetCom{

  my @usage = ("Usage: $0\n

required:
--fastq\tfastq formatted file
--length\tlength of short reads


Will be converted into a fastq file.
	");

	die @usage if ($ARGV[0] eq "");
	GetOptions(\%CMD, "fastq=s", "length=s");

	die("Please specify fastq file\n") unless defined($CMD{fastq});
	die("Please specify length\n") unless defined($CMD{length});
  
	$fastq = $CMD{fastq};
	$length = $CMD{length};


}
	      

