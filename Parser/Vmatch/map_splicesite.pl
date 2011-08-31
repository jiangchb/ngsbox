#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use map_splicesite_hq;

my $mismatches      = shift;
my $read_length     = shift;
my $repeat_mapping  = shift;
my $seedlength      = shift;
my $dir             = shift;
my $suffixtree      = shift;

### Map longer part of spliced read ###
my $half_read_length = $read_length/2;
map_splicesite_hq::map_splicesite($mismatches, $half_read_length, $repeat_mapping, $seedlength, $dir, $suffixtree);


### Extract shorter part of spliced read ###
system( "mkdir $dir/short14 $dir/short15 $dir/short16 $dir/short17 $dir/short18" );
open OUT14, ">$dir/short14/reads_0.fl";
open OUT15, ">$dir/short15/reads_0.fl";
open OUT16, ">$dir/short16/reads_0.fl";
open OUT17, ">$dir/short17/reads_0.fl";
open OUT18, ">$dir/short18/reads_0.fl";

open MAP, "$dir/map.vm" or die "Cannot open $dir/map.vm\n";

while(<MAP>) {
	chomp;
	my @elements = split(/\t/, $_);
	if($elements[5] <= 22) {
		if( ($elements[3] eq "D") && ($elements[6] == 0) ) {
			my $short_end = substr($elements[11], $elements[5], length($elements[11]) );
			if   (length($short_end) == 14 ) {print OUT14 "$elements[0]\t$short_end\n";}
			elsif(length($short_end) == 15 ) {print OUT15 "$elements[0]\t$short_end\n";}
			elsif(length($short_end) == 16 ) {print OUT16 "$elements[0]\t$short_end\n";}
			elsif(length($short_end) == 17 ) {print OUT17 "$elements[0]\t$short_end\n";}
			elsif(length($short_end) == 18 ) {print OUT18 "$elements[0]\t$short_end\n";}
			else { print STDERR "ERROR\n"; }
		}
	}
}

close MAP; close OUT14; close OUT15; close OUT16; close OUT17; close OUT18;


### Map shorter part of spliced read ###
map_splicesite_hq::map_splicesite($mismatches, 14, $repeat_mapping, $seedlength, "$dir/short14", $suffixtree);
map_splicesite_hq::map_splicesite($mismatches, 15, $repeat_mapping, $seedlength, "$dir/short15", $suffixtree);
map_splicesite_hq::map_splicesite($mismatches, 16, $repeat_mapping, $seedlength, "$dir/short16", $suffixtree);
map_splicesite_hq::map_splicesite($mismatches, 17, $repeat_mapping, $seedlength, "$dir/short17", $suffixtree);
map_splicesite_hq::map_splicesite($mismatches, 18, $repeat_mapping, $seedlength, "$dir/short18", $suffixtree);


### Find read combination ###
open MAP, "$dir/map.vm" or die "Cannot open $dir/map.vm\n";
open OUT, ">$dir/halleluja.txt" or die "Cannot open $dir/halleluja.txt\n";

open MAP14, "$dir/short14/map.vm";
open MAP15, "$dir/short15/map.vm";
open MAP16, "$dir/short16/map.vm";
open MAP17, "$dir/short17/map.vm";
open MAP18, "$dir/short18/map.vm";

my $line_14 = <MAP14>;
my $line_15 = <MAP15>;
my $line_16 = <MAP16>;
my $line_17 = <MAP17>;
my $line_18 = <MAP18>;

my @elements_14 = split(/\t/, $line_14);
my @elements_15 = split(/\t/, $line_15);
my @elements_16 = split(/\t/, $line_16);
my @elements_17 = split(/\t/, $line_17);
my @elements_18 = split(/\t/, $line_18);

while(<MAP>) {
	chomp;
	my @elements = split(/\t/, $_);
	if($elements[5] <= 22) {
		if( ($elements[3] eq "D") && ($elements[6] == 0) ) {
			my $short_end = substr($elements[11], $elements[5], length($elements[11]) );
                        if(length($short_end) == 14 ) {
                                while($elements_14[0] eq $elements[0]) {
                                        if(     ($elements_14[3] eq "D") &&
                                                ($elements_14[1] == $elements[1]) &&
                                                (abs($elements_14[2] - $elements[2]) < 20000)
                                        ) {
                                                print OUT "$elements[0]\t$elements[1]\t$elements[2]\t" .
                                                        "$elements[7]\t$elements_14[1]\t$elements_14[2]\t" .
                                                        "$elements_14[7]\t$elements[8]\n";
                                        }
                                        $line_14 = <MAP14>;
                                        @elements_14 = split(/\t/, $line_14);
                                }
                        }
			elsif(length($short_end) == 15 ) {
				while($elements_15[0] eq $elements[0]) {
					if( 	($elements_15[3] eq "D") && 
						($elements_15[1] == $elements[1]) &&
						(abs($elements_15[2] - $elements[2]) < 20000)
					) {
						print OUT "$elements[0]\t$elements[1]\t$elements[2]\t" .
							"$elements[7]\t$elements_15[1]\t$elements_15[2]\t" .
							"$elements_15[7]\t$elements[8]\n";
					}
					$line_15 = <MAP15>;
					@elements_15 = split(/\t/, $line_15);
				}
			}
			elsif(length($short_end) == 16 ) {
				while($elements_16[0] eq $elements[0]) {
					if(     ($elements_16[3] eq "D") &&
						($elements_16[1] == $elements[1]) &&
						(abs($elements_16[2] - $elements[2]) < 20000)
					){
						print OUT "$elements[0]\t$elements[1]\t$elements[2]\t" .
							"$elements[7]\t$elements_16[1]\t$elements_16[2]\t" .
							"$elements_16[7]\t$elements[8]\n";
					}
					$line_16 = <MAP16>;
					@elements_16 = split(/\t/, $line_16);
				}
			}
			elsif(length($short_end) == 17 ) {
				while($elements_17[0] eq $elements[0]) {
					if(     ($elements_17[3] eq "D") &&
						($elements_17[1] == $elements[1]) &&
						(abs($elements_17[2] - $elements[2]) < 20000)
					){
						print OUT "$elements[0]\t$elements[1]\t$elements[2]\t" .
							"$elements[7]\t$elements_17[1]\t$elements_17[2]\t" .
							"$elements_17[7]\t$elements[8]\n";
					}
					$line_17 = <MAP17>;
					@elements_17 = split(/\t/, $line_17);
				}
			}
			elsif(length($short_end) == 18 ) {
				while($elements_18[0] eq $elements[0]) {
					if(     ($elements_18[3] eq "D") &&
						($elements_18[1] == $elements[1]) &&
						(abs($elements_18[2] - $elements[2]) < 20000)
					){
						print OUT "$elements[0]\t$elements[1]\t$elements[2]\t" .
							"$elements[7]\t$elements_18[1]\t$elements_18[2]\t" .
							"$elements_18[7]\t$elements[8]\n";
					}
					$line_18 = <MAP18>;
					@elements_18 = split(/\t/, $line_18);
				}
			}
		}
	}
}

exit(0);
