#! /usr/bin/perl

use strict;

my $usage = "$0 coveragefile genotypingfile markerfile coveragefactor\n";
my $coverage = shift or die $usage; #"marker_coverage.txt";
my $genotyping = shift or die $usage; #"genotyping.txt";
my $markerfile = shift or die $usage; #"Ler-1 scheiss.txt";
my $covfactor = shift or die $usage; # Multipe with col-0 cov;

my %COV = ();
my %REF = ();
my %MARKER = ();

open COVERAGE, $coverage or die "Cannot open file\n";
while (my $line = <COVERAGE>) {
	my @a = split " ", $line;
	$COV{$a[0]."#".$a[1]} = $a[2];
	$REF{$a[0]."#".$a[1]} = $a[3];
}
close COVERAGE;

open MARKER, $markerfile or die "Cannot open file\n";
while (my $line = <MARKER>) {
        my @a = split " ", $line;
        $MARKER{$a[1]."#".$a[2]} = $a[4];
}
close MARKER;



open FILE, $genotyping or die "Cannot open file\n";
while (my $line = <FILE>) {
        my @a = split " ", $line;

	if ($REF{$a[0]."#".$a[1]} ne "N" and defined($COV{$a[0]."#".$a[1]})) {
		my $parent_count = $a[2] + $a[3];
		my $mutant_bg_count = $a[2];

		my $mutant_allele = 0;
		my $second_allele = 0;
		
		for (my $i = 0; $i < ($COV{$a[0]."#".$a[1]} * $covfactor); $i++) {
			my $draw = int(rand($parent_count))+1;

			# Introduce sequencing errors:
			my $seqerr = int(rand(1000))+1;
			if ($seqerr<=3) {
				if ($draw <= $mutant_bg_count) {
					$draw = $mutant_bg_count+1;
				}
				else {
					$draw = $mutant_bg_count-1;
				}
			}

			if ($draw <= $mutant_bg_count) {
				$mutant_allele++;
			}
			else {
				$second_allele++;
			}
		}

		########################################################################
		# Print out file in consensus summary format
		########################################################################
		print $a[0], "\t", $a[1], "\t-\t", $COV{$a[0]."#".$a[1]}; 

		if ($REF{$a[0]."#".$a[1]} eq "A") {
			print "\t", $mutant_allele; 
		}
		elsif ($MARKER{$a[0]."#".$a[1]} eq "A")  {
			print "\t", $second_allele;
		}
		else {
			print "\t0";
		}

		if ($REF{$a[0]."#".$a[1]} eq "C") {
                        print "\t", $mutant_allele; 
                }
		elsif ($MARKER{$a[0]."#".$a[1]} eq "C")  {
                        print "\t", $second_allele;
                }
                else {
                        print "\t0";
                }


		if ($REF{$a[0]."#".$a[1]} eq "G") {
                        print "\t", $mutant_allele;
                }
		elsif ($MARKER{$a[0]."#".$a[1]} eq "G")  {
                        print "\t", $second_allele;
                }
                else {
                        print "\t0";
                }

                if ($REF{$a[0]."#".$a[1]} eq "T") {
                        print "\t",$mutant_allele;
                }
		elsif ($MARKER{$a[0]."#".$a[1]} eq "T")  {
                        print "\t", $second_allele;
                }
                else {
                        print "\t0";
                }

		for (my $i = 0; $i < 36; $i++) {
			print "\t0";
		}

		print "\t", $REF{$a[0]."#".$a[1]};

		print "\n";
		########################################################################

	}
}
	




