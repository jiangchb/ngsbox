#! /usr/bin/perl
use strict;
use Getopt::Long;

my %CMD;
my $file;
my $out;
my $idx;

GetCom();

open FILE, $file or die "Need a fasta file\n";
open OUT, ">".$out or die "Need a fasta file\n";
open IDX, ">".$idx or die "Need a fasta file\n";

my $count = 1;
while (<FILE>) {
	chomp();
	if (substr($_, 0, 1) eq ">") {
		my @a = split " ", $_;
		print OUT ">$count\n";
		#print IDX $count, "\t", substr($a[0], 1, length($a[0])-1), "\n";
		print IDX $count, "\t", substr($_, 1, length($_)-1), "\n";
		$count++;
	} else {
		print OUT $_, "\n";
	}
}


sub GetCom{

        my @usage = ("\nUsage: $0 --file=file --out=filename --idx=filename
                --file\treference genome
		--out\treference genome with index in fasta headers (used for mapping)
		--idx\t1 to 1 look-up for reference file (original fasta entry to idx fasta entry)

                description:

                \n\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "out=s", "idx=s");

        die("Please specify file\n")    if not defined $CMD{file};
	die("Please specify out\n")    if not defined $CMD{out};
	die("Please specify idx\n")    if not defined $CMD{idx};

        $file = $CMD{file};
	$out = $CMD{out};
	$idx = $CMD{idx};

}



