#!/usr/bin/perl

use strict;
use warnings;


### User params
my $usage = "\n$0 ecotype assembly reference chromosome maxDivBases outfolder\n\n";

my $ecotype   = shift or die $usage;
my $assembly  = shift or die $usage;
my $reference = shift or die $usage;
my $chr       = shift or die $usage;
my $maxDivBp  = shift or die $usage;
my $outfolder = shift or die $usage;


### Init other variables
my $chr_int = substr($chr, 3);


### Create folders
mkdir "$outfolder" or die;
mkdir "$outfolder/Alignments" or die;


### Run NUCMER
system("nucmer --mum -b $maxDivBp -g 90 -l 35 -c 80 --prefix=$outfolder/wgaMum $reference $assembly");

### Run delta filter
system("delta-filter -q $outfolder/wgaMum.delta > $outfolder/wgaMum.filter");


### Parse reference file
open REFERENCE, $reference or die "Cannot open $reference\n";

my $id = "";
my $REF = "";

while(<REFERENCE>) {
        chomp;
	
	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		if($id eq $chr) {
			$REF .= $_;
		}
	}
}


### Parse assembly file
open ASSEMBLY, $assembly or die "Cannot open $assembly\n";

$id = "";
my %CTG = ();

while(<ASSEMBLY>) {
	chomp;

	if (substr($_, 0, 1) eq ">") {
		my $header = substr($_, 1);
		my @e = split(" ", $header);
		$id = $e[0];
	}
	else {
		$CTG{$id} .= $_;
	}
}
close ASSEMBLY;



### Open result files and hashes
open INVFILE, ">$outfolder/inversion.txt" or die  "Cannot open $outfolder/inversion.txt\n";

### Create alignments for each scaffold ID and parse HDRs and Inversions
foreach $id (keys %CTG) {

	### Create and open alignment for one scaffold
	system("show-aligns -r $outfolder/wgaMum.filter $chr $id > $outfolder/Alignments/alignments_$id.txt");
	open ALIGNMENT, "$outfolder/Alignments/alignments_$id.txt" or die "Cannot open $outfolder/Alignments/alignments_$id.txt\n";


	### Init variables
	my $state = 0;
	my $junk_line = "";
	my $refseq = "";
	my $scaffseq = "";
	my $HSPcounter = 0;
	my %HSP = ();
	my %INV = ();


	### Parse alignment file of one scaffold
	while(<ALIGNMENT>) {
		chomp;

		### New HSP found ----------------------------------------------------------------------
		if( $_ =~ /BEGIN/ ) {
			my $last_HSP_counter = $HSPcounter;
			$HSPcounter++;
			my @e = split(" ", $_);
			my %current_HSP = ( 'refbeg' => $e[5], 'refend' => $e[7], 'scaffbeg' => $e[10], 'scaffend' => $e[12] );
			$HSP{$HSPcounter} = \%current_HSP;


			# Check if data for inversion of last HSP has to be completed
			if(exists $INV{$last_HSP_counter}) {
				$INV{$last_HSP_counter}{next_refbeg} = $e[5];
				$INV{$last_HSP_counter}{next_scaffbeg} = $e[10];

				if(	( $INV{$last_HSP_counter}{next_refbeg} > $INV{$last_HSP_counter}{refend} - 20) &&
					( $INV{$last_HSP_counter}{next_refbeg} - $INV{$last_HSP_counter}{refend} <= 1000) &&
					( $INV{$last_HSP_counter}{next_scaffbeg} > $INV{$last_HSP_counter}{scaffbeg} - 20) &&
					( $INV{$last_HSP_counter}{next_scaffbeg} - $INV{$last_HSP_counter}{scaffbeg} <= 1000 )
				){
					my $len = $INV{$last_HSP_counter}{scaffbeg} - $INV{$last_HSP_counter}{scaffend} + 1;
					
					print INVFILE "$ecotype\tINV\t$chr_int\t$id\t$last_HSP_counter\t$len\t" .
						$INV{$last_HSP_counter}{last_refend} ."\t". $INV{$last_HSP_counter}{refbeg} ."\t". 
						$INV{$last_HSP_counter}{refend} ."\t". $INV{$last_HSP_counter}{next_refbeg} ."\t".
						$INV{$last_HSP_counter}{last_scaffend} ."\t". $INV{$last_HSP_counter}{scaffbeg} ."\t". 
						$INV{$last_HSP_counter}{scaffend} ."\t". $INV{$last_HSP_counter}{next_scaffbeg} ."\n";
				}
			}

			# Check if scaffold is reversed in HSP alignment
			if($e[12] < $e[10]) {

				# Check minimum size of inversion and upstream/downstream anchors
				if( ($e[7] - $e[5] > 50) && ($e[10] - $e[12] > 50) ) { #&& ( $e[12] > 100 ) ) {

					# Check if reference and scaffold positions are consistent
					if( 	(exists $HSP{$HSPcounter - 1} ) &&
						($HSP{$HSPcounter - 1}{refend} < $e[5] + 20) && 
						($e[5] - $HSP{$HSPcounter - 1}{refend} <= 1000) &&
						($HSP{$HSPcounter - 1}{scaffend} < $e[12] + 20) &&
						($e[12] - $HSP{$HSPcounter - 1}{scaffend} <= 1000)
					) { 

						my %current_INV = ( 	last_refend => $HSP{$HSPcounter - 1}{refend},
									refbeg => $e[5],
									refend => $e[7],
									next_refbeg => -1,
									last_scaffend => $HSP{$HSPcounter - 1}{scaffend},
									scaffbeg => $e[10],
									scaffend => $e[12],
									next_scaffbeg => -1
								);

						$INV{$HSPcounter} = \%current_INV;
					}
				}
			}
			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1;
		}


		### HSP finished, process information --------------------------------------------------
		elsif(  $_ =~ /END/ ) {

			### Reset
			$refseq = "";
			$scaffseq = "";
			$state = 0;

			### compare sequence of inverted sequence
			# TODO
		}


		### Correct for stupid extra empty line after HSP-END ----------------------------------
		elsif( $_ eq "" ) {
			# Do nothing
		}


		### Reference line of alignment --------------------------------------------------------
		elsif( $state == 1 ) {
			my ($refpos, $seq) = split(" ", $_);
			$refseq .= $seq;
			$state = 2;
		}


		###Scaffold line of alignment ----------------------------------------------------------
		elsif( $state == 2 ) { 
			my ($scaffpos, $seq) = split(" ", $_);
			$scaffseq .= $seq;
			$junk_line = <ALIGNMENT>;
			$junk_line = <ALIGNMENT>;
			$state = 1; 
		}
	}

	close ALIGNMENT;
}

close INVFILE;

print STDERR "Finished Whole Genome Alignment\n";
exit(0);
