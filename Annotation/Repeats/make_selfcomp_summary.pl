#!/usr/bin/perl
###############################################################
#Author 	Korbinian, Stephan Ossowski
#Date 		07/03/07
#Version	0.9
#Input		Mapping of self comparison, repeat summary
###############################################################

use strict;
use warnings;
use File::Temp;
use Cwd;


### Parse input
my $file_frags    = shift;

if( (not defined $file_frags)  ) {
	print "Usage example: perl make_selfcomp_summary.pl map.list > consensus.out"; 
}

open (FRAGS, $file_frags) or die "Cannot open input file $file_frags\n";

### Global variables
my $current_chr        = 1;
my $current_pos        = 1;

### Input arrays
my @matrix             = ();	# short read sequence matrix
my @mismatches         = ();	# number of mismatches per read
my @hits_array         = ();	# repetitiveness of reads
my @pos                = ();	# positions in reads
my @read_length        = ();	# read length (variance due to end gap)
my @orientations       = ();	# orientation array
my @errors_per_read    = ();	# estimated number of errors per read
my %insertions         = ();	# hash with potential insert position

### Mapping features
# General counts
my %count_base                 = (A => 0, C => 0, G => 0, T => 0, '-' => 0, N=> 0);
my @mismatch_type              = (0, 0, 0, 0, 0);	# MM0, MM1, MM2, MM3, MM4
my $coverage                   = 0;
my $sum_hits                   = 0;
my $sum_mismatches             = 0;
my $call                       = "";

###########################################################################################################################
### Scan input file
while( <FRAGS> ) {
	chomp($_);
	my @members = split(/\t/, $_);
	my $chr            = $members[0];
	my $position       = $members[1];
	my $s              = $members[2];
	my $solexa_id      = $members[3];
	my $orientation    = $members[4];
	my $mismatch       = $members[5];
	my $hits           = $members[6];
	my $used_length    = $members[7];
	my $offset         = $members[8];

	### New position reached, sum up features of last position ########################################################
	if( $position != $current_pos ) {
		my $next = $position;
		if($chr != $current_chr) { $next = $current_pos + length($matrix[0]);}

		for(my $i = $current_pos; $i < $next; $i++) {
			my @analysis_base     = ();
			my @analysis_pos      = ();
			my @analysis_length   = ();
			my @analysis_orientation = ();
			#print DEBUG "$position\t$current_pos\t$next\t$i\t" . scalar(@matrix) . "\n";

			### Leave loop if nothing is in matrix
			if ( scalar(@matrix) == 0 ) {
				#print DEBUG "$position\t$current_pos\t$next\t$i\n";
				#print 	"$current_chr\t$i" .
				#	"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" .
				#	"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" .
				#	"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t" .
				#	"\t\t";
				next;
			}

			### Run through all bases (reads) overlapping the current position
			for(my $j = scalar(@matrix)-1; $j >= 0; $j--) {
				my $base_call = substr($matrix[$j], 0, 1);
				unshift @analysis_base, $base_call;
				unshift @analysis_pos, $pos[$j];
				unshift @analysis_length, $read_length[$j];
				unshift @analysis_orientation, $orientations[$j];

				### Update general counts
				$coverage++;
				$sum_hits += $hits_array[$j];
				$sum_mismatches += $mismatches[$j];
				$mismatch_type[$mismatches[$j]]++;
				$count_base{$base_call}++;
				
				### Delete first character of string
				substr($matrix[$j], 0, 1) = "";

				### Increment genome position
				$pos[$j]++;
				
				### Delete empty strings
				if($matrix[$j] eq "") {
					splice(@matrix, $j, 1);
					splice(@mismatches, $j, 1);
					splice(@hits_array, $j, 1);
					splice(@pos, $j, 1);
					splice(@read_length, $j, 1);
					splice(@orientations, $j, 1);
					splice(@errors_per_read, $j, 1);
				}
			}


			### Determine called base
			my @sorted_call   = reverse sort { $count_base{$a} <=> $count_base{$b} } keys %count_base;
			$call             = $sorted_call[0];
				
			### Calculate average mismatches and average hits (repetitiveness)
			my $average_mismatches = $sum_mismatches / $coverage;
			my $average_hits = $sum_hits / $coverage;
			
			### Print summary file (consensus stage 1)
			print   "$current_chr\t$i\t";

			print 	"$call\t$coverage\t" .
				"$count_base{A}\t$count_base{C}\t$count_base{G}\t" .
				"$count_base{T}\t$count_base{'-'}\t$count_base{N}\t" .
				sprintf("%.2f",$average_hits) . "\t".
				sprintf("%.2f",$average_mismatches) . "\t".
				"$mismatch_type[0]\t$mismatch_type[1]\t$mismatch_type[2]\t" .
				"$mismatch_type[3]\t$mismatch_type[4]\n";

			### Reset variables
			#General counts
			%count_base                 = (A => 0, C => 0, G => 0, T => 0, '-' => 0, N => 0);
			@mismatch_type              = (0, 0, 0, 0, 0);
			$coverage                   = 0;
			$sum_hits                   = 0;
			$sum_mismatches             = 0;

		}
	}




	### Do for each sequence #######################################################################################

        my $posinref       = $position;
        my $seq            = "";
        my $i              = 0;

	### Find insertions in read, create seq without brackets and inserts
        while($i < length($s)) {
                my $char = substr($s, $i, 1);

                # Mismatch 
                if($char eq "[") {
                        $i++;
                        $char = substr($s, $i, 1);
                        $i++;
                        $char = substr($s, $i, 1);
                        $posinref++;
                        $seq .= $char;
                        $i += 2;
                }
                else {
                        $posinref++;
                        $seq .= $char;
                        $i++;
                }
        }


	$seq =~ s/I//g;

        ### Add new sequence to the beginning of the matrix
        unshift(@matrix, $seq);
        unshift(@mismatches, $mismatch);
	unshift(@errors_per_read, 0);
        unshift(@hits_array, $hits);
        unshift(@pos, 0);
        unshift(@read_length, length($seq));
	unshift(@orientations, $orientation);

	### update genome position
	$current_pos = $position;
	$current_chr = $chr;
}




exit(0);

