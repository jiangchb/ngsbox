#!/usr/bin/perl

# --------------------------------------------------------------------------
# Build a genome matrix of hundreds of transcriptomes
# Written by Stephan Ossowski
# --------------------------------------------------------------------------

use strict;
use warnings;

use Getopt::Long;
use FindBin;

### Command line parameters -------------------------------------------------------------------
my $num_samples = "";
my $refseq      = "";
my $folder_list = "";

my %CMD;
GetCom();

my %coverage = ();


### Read reference sequence
my $chr = "";
my $seq = "";
open REFSEQ, $refseq or die "Cannot open $refseq\n";
while( <REFSEQ> ) {
        chomp;
        if($_ =~ /^>/ ) {
                if($seq ne "") {
                        for(my $i = 1; $i <= length($seq); $i++) {
				my @tmp = ();
				for(my $x = 0; $x < $num_samples; $x++) {
					push( @tmp, 0 );
				}
				$coverage{"$chr#$i"} = \@tmp;
                        }
                }

                ($chr) = split(" ", $_);
                $chr =~ s/>//g;
                $seq = "";
        }
        else {
                $seq .= $_;
        }
}
if($seq ne "") {
        for(my $i = 1; $i <= length($seq); $i++) {
                my @tmp = ();
		for(my $x = 0; $x < $num_samples; $x++) {
			push( @tmp, 0 );
		}
		$coverage{"$chr#$i"} = \@tmp;
        }
}
$seq = "";
close REFSEQ;


### Prepare database tables and ecotype list
my $create_table = "CREATE table coverage_matrix( chromosome tinyint(1), position int(9)";
my $eco_list = "";
my $counter = 0;

### Read all result files from project folders
open FOLDER, $folder_list or die "Cannot open $folder_list\n\n";
while( <FOLDER> ) {
	chomp;
	my $folder_name = $_;

	### Open consensus file
	open CONS, "$folder_name/consensus_summary.txt" or next;


	### Prepare database tables and ecotype list
	my @folder_structure = split("\/", $folder_name);
	my $eco = $folder_structure[0];
	$eco_list .= ($counter + 3) . "\t$eco\n";
	$create_table .= ", " . $eco . "_all_coverage int(6)";


	### Parse files and store coverage (all and nonrep)
	while( <CONS> ) {
		chomp;
		my ($chr, $pos, $ref, $cov) = split("\t", $_);
		#$coverage{$a[0]."#".$a[1]}[$counter] = $a[3];
		$coverage{"$chr#$pos"}[$counter] = $cov;
	}

	### Finish processing file
	close CONS;
	$counter++;
}


### Print coverage frequency matrix
open OUT, ">coverage_matrix.txt" or die "Cannot open output file\n";
foreach my $chrpos (keys %coverage) {
	my ($chr, $pos) = split("#", $chrpos);
	print OUT "$chr\t$pos";
	for(my $i = 0; $i < $counter; $i++) {
		print OUT "\t" . $coverage{$chrpos}[$i];
	}
	print OUT "\n";
}


### Print database create statement
$create_table .= ", PRIMARY KEY(chromosome, position));";
open CREATE, ">coverage_matrix_create_table.txt" or die "Cannot open create file\n";
print CREATE "$create_table";
close CREATE;


### Print list of ecotypes
open LIST, ">coverage_matrix_ecotype_list.txt" or die "Cannot open ectype file\n";
print LIST $eco_list;
close LIST;


### Finish
exit(0);



### Read command line parameters --------------------------------------------------------------
sub GetCom {
  my @usage = ("\nUsage: $0

--samples    INT        Number of samples t be processed
--refseq     STRING     Reference sequence file
--folder     STRING     Text file with folder list of resequencing projects
\n");


	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "samples=s", "refseq=s", "folder=s");

	die("Please specify number of samples\n") unless defined($CMD{samples});
        die("Please specify ref-seq file\n") unless defined($CMD{refseq});
	die("Please specify folder file\n") unless defined($CMD{folder});

	$num_samples = $CMD{samples};
	$refseq      = $CMD{refseq};
	$folder_list = $CMD{folder};
}
