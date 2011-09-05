#! /usr/bin/perl
use strict;
use warnings;

###### 
# NGSbox - bioinformatics analysis tools for next generation sequencing data
#
# Copyright 2007-2011 Stephan Ossowski, Korbinian Schneeberger
# 
# NGSbox is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.
#
# NGSbox is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# Please find the GNU General Public License at <http://www.gnu.org/licenses/>.
#
#  -------------------------------------------------------------------------
#
#  Module: Parser::ML::GraphMapping::graphmap_correction.pl
#  Purpose:
#  In:
#  Out:
#



# The transformation of indels has errors
# This scripts tries to correct them.

my $usage = "$0 map.list\n";
my $file = shift or die $usage;
open FILE, $file or die $usage;

while (my $line = <FILE>) {
	my @a = split " ", $line;
	my $eco = $a[0];
	if ($eco eq "Col-0") {
		if ($a[4] eq $a[5]) {
			print $line;
		}
		else {
			print $a[0], "\t", $a[1], "\t", $a[2], "\t", $a[3], "\t", $a[5], "\t", $a[5];
			for (my $i = 6; $i < @a+0; $i++) {
				print "\t", $a[$i];
			}
			print "\n";
		}
	}
	else  {
		######################################################
		# Correct for Gap MM Error:
		my $num_gap_mm_error_corr = 0;
		my $new = "";
		if ($a[4] =~ m/\[.*\]\(.*\-/) { # only pre-selection
			for (my $i = 0; $i < length($a[4]); $i++) {
				if (substr($a[4], $i, 1) eq "[" && substr($a[4], $i+3, 1) eq "]" && substr($a[4], $i+4, 1) eq "(" && substr($a[4], $i+7, 1) eq ")" && substr($a[4], $i+6, 1) eq "-" && substr($a[4], $i+2, 1) eq substr($a[4], $i+5, 1)) {
					$new .= "[".substr($a[4], $i+1, 1)."-]".substr($a[4], $i+2, 1);
					$i+=7;
					$num_gap_mm_error_corr++;
				}
				elsif (substr($a[4], $i, 1) eq "[" && substr($a[4], $i+3, 1) eq "]" && substr($a[4], $i+4, 1) eq "(" && substr($a[4], $i+7, 1) eq ")" && substr($a[4], $i+5, 1) eq "-" && substr($a[4], $i+1, 1) eq substr($a[4], $i+6, 1)) {
					$new .= "[-".substr($a[4], $i+2, 1)."]".substr($a[4], $i+1, 1);
					$i+=7;
					$num_gap_mm_error_corr++;
				}
				else {
                                        $new .= substr($a[4], $i, 1);
                                }
			}
		}

		$a[8] -= $num_gap_mm_error_corr;
		if ($new ne "") {
			$a[4] = $new;
		}

		######################################################
		# Correct for most left gap
		$num_gap_mm_error_corr = 0;
		if (($a[4] =~ m/\[-.\]/) or ($a[4] =~ m/\(-.\)/)) { # only pre-selection
			my $changed = 0;
			do {
				$new = "";	
				$changed = 0;
				for (my $i = 0; $i < length($a[4]) && $changed == 0; $i++) {
					# Strore i
					my $j = $i;
					# detect full length of insertion
					my $ins = "";
					while ((substr($a[4], $i, 1) eq "[" or substr($a[4], $i, 1) eq "(") and substr($a[4], $i+1, 1) eq "-") {
						$ins .= substr($a[4], $i+2, 1);
						$i+=4;
					}
					# Check if length(ins) bases no indel, but MM occurs
					my $mm_detected = 0;
					my $ref = "";
					my $post_ins;
					if ($ins ne "") {
						for (my $k = 0; $k < length($ins) && $mm_detected >= 0; $k++) {
							if (substr($a[4], $i, 1) eq "[" or substr($a[4], $i, 1) eq "(") {
								if (substr($a[4], $i+1, 1) eq "-" or substr($a[4], $i+2, 1) eq "-") {
									$mm_detected = -1;
								}
								else {
									$mm_detected++;
								}
								$ref .= substr($a[4], $i+1, 1);
								$post_ins .= substr($a[4], $i+2, 1);
								$i += 4;
							}
							else {
								$ref .= substr($a[4], $i, 1);
								$post_ins .= substr($a[4], $i, 1);
								$i++;
							}
						}
					}
					# Reset alignment
					if ($mm_detected == 1 and $ref eq $ins) {
						$num_gap_mm_error_corr += $mm_detected;
						$new = substr($a[4], 0, $j) . $ins;
						for (my $r = 0; $r < length($post_ins); $r++) {
							$new .= "[-".substr($post_ins, $r, 1)."]";
						}
						$new .= substr($a[4], $i, length($a[4])-$i); 
						$changed = 1;
					}
					
					# Reset i
					$i = $j;
				}

                		if ($new ne "") {
		                        $a[4] = $new;
					# Set number of mm
					my $count_mm = 0;
					for (my $s = 0; $s<length($a[4]); $s++) {
						if (substr($a[4], $s, 1) eq "[" or substr($a[4], $s, 1) eq "(") {
							$s+=3;
							$count_mm++;
						}
					}
					$a[8] = $count_mm;
				}


			} while ($changed == 1) ;
		}

			
		######################################################
		# Write corrected mapping:
		for (my $i = 0; $i < @a; $i++) {
                	print $a[$i];
			print "\t" if $i != @a-1;
                }
		print "\n";
	}
}

