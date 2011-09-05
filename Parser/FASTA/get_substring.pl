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
#  Module: Parser::FASTA::get_substring.pl
#  Purpose:
#  In:
#  Out:
#



my $usage = "$0 file newFastaID id begin end\nIf begin is larger than end, sequence will be reported rev comp from end to begin\n";

my $file   = shift or die $usage;
my $new_id = shift or die $usage;
my $id     = shift or die $usage;
my $begin  = shift or die $usage;
my $end    = shift or die $usage;

open FILE, $file or die "Cannot open $file\n";

my $seq = "";
my $flag = 0;

while(my $line = <FILE>) {
	chomp($line);
	my @a = split " ", $line;
	if (substr($a[0], 0, 1) eq ">") {
		if ($flag == 1) {
			if ($begin <= $end) {
				print ">".$new_id." ".$id."_".$begin."_".$end."\n";
				print substr($seq, $begin-1, $end-$begin+1), "\n";
			}
			else {
				print ">".$new_id." ".$id."_".$begin."_".$end."_revcomp\n";
				my $seq = substr($seq, $end-1, $begin-$end+1);
				$seq = revcomp($seq);
				print $seq, "\n";
			}
			exit(0);
		}

		if ($id eq substr($a[0], 1, length($a[0])-1)) {
			$flag = 1;
		}
		else {
			$flag = 0;
		}

                $seq = "";
	}
	else {
		$seq .= $line;
	}
}

if ($flag == 1) {
	if ($begin <= $end) {
 	       print ">".$new_id." ".$id."_".$begin."_".$end."\n";
        	print substr($seq, $begin-1, $end-$begin+1), "\n";
         }
                        else {
                                print ">".$new_id." ".$id."_".$begin."_".$end."_revcomp\n";
                                my $seq = substr($seq, $begin-1, $end-$begin+1);
                                $seq = revcomp($seq);
                                print $seq, "\n";
                        }
}

exit(0);


sub revcomp {
        my ($seq) = @_;

        my $newseq = reverse $seq;
        $newseq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;

        return $newseq;
}

