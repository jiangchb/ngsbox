#! /usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;

my %CMD;

my $hist = 0;
my $verbose = 0;

GetCom();

my $a = 0;
my $c = 0;
my $g = 0;
my $t = 0;
my $num = 0;

while (<FILE>) {
	my @a = split " ", $_;
	for (my $i=0; $i<length($a[1]); $i++) {
		my $char = substr($a[1], $i, 1);
		$num++;
		if ($char eq "A" or $char eq "a") {
			$a++;
		}
		if ($char eq "C" or $char eq "c") {
        	        $c++;
	        }
        	if ($char eq "G" or $char eq "g") {
                	$g++;
	        }
        	if ($char eq "T" or $char eq "t") {
                	$t++;
	        }
	}
	print STDERR $num if $num%500000 == 0 and $verbose == 1;
}

print "Base type count of ".$CMD{file}.":\n";
print "A: ", $a, "\t", $a/$num,"\n";
print "C: ", $c, "\t", $c/$num,"\n";
print "G: ", $g, "\t", $g/$num,"\n";
print "T: ", $t, "\t", $t/$num,"\n";
print "\n";
print "GC: ", $c+$g, "\t", ($c+$g)/$num, "\n";
print "AT: ", $a+$t, "\t", ($a+$t)/$num, "\n";

sub GetCom{

  my @usage = ("Usage: base_type_dist_from_fl.pl --file=<fl file> -v\n
required:
--file\t\tDefine the fl file for GC content calculation
optional:
-v\t\tReport every 10.000th line
\n");

        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "v");

        die("Please specify a call string\n") unless $CMD{file};
	
        if (defined($CMD{v})) {
                $verbose = 1;
        }


	open FILE, $CMD{file} or die "Cannot open file\n";
}




