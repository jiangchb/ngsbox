#! /usr/bin/perl

my $usage = "perl parse_agf.pl afg_file input_fasta\n";
open FILE, $ARGV[0] or die $usage;
open FASTA, $ARGV[1] or die $usage;

my $count = 1;
my %READ_ID = ();
while (my $line = <FASTA>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		$READ_ID{$count} = substr($line, 1, length($line)-1);
		$count++;
	}
}

print STDERR "Parsed fasta\n";

my $id;
my $first_flag = 1;
CONTIGS: while (my $line = <FILE>) {
	chomp($line);
	if ($line eq "{CTG") {
		$id = <FILE>;
		chomp($id);
		$id =~ s/iid://g;
		if ($first_flag == 0) {
			print "\n";
		}
		else {
			$first_flag = 0;
		}
		print $id;
		<FILE>; <FILE>;
		my $seq = "";
		SEQ: while (my $l = <FILE>) {
			chomp($l);
			if (substr($l, 0, 1) eq ".") {
				print "\t", length($seq), "\t", $seq;
				last SEQ;	
			}
			$seq .= $l;	
		}
	}
	elsif ($line eq "{TLE") {
		$tile = <FILE>;
		chomp($tile);
		$tile =~ s/src://g;
		print "\t", $READ_ID{$tile};
	}
	elsif ($line eq "{RED") {
		last CONTIGS;
	}
}



