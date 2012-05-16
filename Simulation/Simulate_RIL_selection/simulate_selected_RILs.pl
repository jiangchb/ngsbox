#! /usr/bin/perl
use strict;

use Math::Random qw(random_normal);

my $usage = "$0 run_id\n";
my $RUN_ID = shift or die $usage;

my $QTL_NUM = 1;
my $HERITABILITY = 0.5;
my $EXTREME_PERCENTILE = 0.1;
my $POPULATION_SIZE = 10000;
#my $POPULATION_FILE = "Simulation_input/genotyped_f5_5000_merged";
my $POPULATION_FILE = "Simulation_input/genotyped_f5_33k_merged";
my $MARKER_FILE = "Simulation_input/Ler-1.331000.txt";
my $EFFECT_FILE = "Simulation_input/".$QTL_NUM."QTL_geodist.txt";
my $OUTPUT_FOLDER = "Simulated_pools";

###########################################################################
## write log file

open LOG, "> $OUTPUT_FOLDER/Sim.$RUN_ID.log" or die "cannot open log file";
print LOG"#run_id\tqtl_num\theritability\textreme_percentile\tpop_size\tpop_file\tmarker_file\teffect_file\toutput_folder\n";
print LOG "$RUN_ID\t$QTL_NUM\t$HERITABILITY\t$EXTREME_PERCENTILE\t$POPULATION_SIZE\t$POPULATION_FILE\t$MARKER_FILE\t$EFFECT_FILE\t$OUTPUT_FOLDER\n";

############################################################################
# Read in effects of the QTLs
my %QTL_EFF = ();
open FILE, $EFFECT_FILE or die "cannot find the respective effects file\n";
while (<FILE>) {
        my @a = split, " ";
        my $qtl_id = $a[0];
        $qtl_id =~ s/qtl//g;
        if ($a[1] == 1) {
                $QTL_EFF{$qtl_id}{"C"} = $a[2];
        }
        else {
                $QTL_EFF{$qtl_id}{"L"} = $a[2];
        }
        $QTL_EFF{$a[0]}{"H"} = 0;
}
close FILE;

print STDERR  "Got effects\n";


###########################################################################
# Read in markers and select some of them as QTL
my %QTL_LOC = ();

open QTL, "> $OUTPUT_FOLDER/Sim.$RUN_ID.qtl.txt" or die "cannot open log file";

open FILE, $MARKER_FILE or die "Cannot find file $MARKER_FILE\n";
my $c = 0;
my $qtl_id = 1;
my $num_marker_remaining = 331729;
my $num_QTL_remaining = $QTL_NUM;
print QTL "#qtl_id\tmarker_id\tchr\tpos\teffect_p1\teffect_p2\n";
while (<FILE>) {
        my @a = split, " ";

	# Randomize if this marker is one of the QTL 
	my $rand = rand();
        if ($rand <= ($num_QTL_remaining / $num_marker_remaining)) {
                $num_QTL_remaining--;
                $QTL_LOC{$qtl_id} = $c;
		print QTL $qtl_id, "\t", $c, "\t", $a[1], "\t", $a[2], "\t", $QTL_EFF{$qtl_id}{"C"}, "\t",  $QTL_EFF{$qtl_id}{"L"}, "\n";
		$qtl_id++;
        }

	$num_marker_remaining--;
	$c++;
}
close FILE;

print STDERR  "Got markers and selected QTLs on it. ".($qtl_id-1)." selected.\n";


###########################################################################
# Read in Genotypes and phenotype them based on the genotypes

my %ID2GENOPHENO = ();

my $SUM_PHENOTYPES = 0;
my $SQUARED_SUM_PHENOTYPES = 0;
my $IND_COUNT = 0;

my $num_remaining = `wc -l $POPULATION_FILE`;
chomp($num_remaining);
$num_remaining--; # header line

print STDERR "$POPULATION_FILE has $num_remaining individuals. $POPULATION_SIZE will be phenotyped.\n";
if ($num_remaining < $POPULATION_SIZE) {
	die("Too little genotypes in file..\n");
}

my $num_ind_remaining = $POPULATION_SIZE;

open FILE, $POPULATION_FILE or die "cannot open $POPULATION_FILE\n";
my $header = <FILE>;
while (my $l = <FILE>) {
	$IND_COUNT++;
	print STDERR $IND_COUNT, "\n" if ($IND_COUNT%100 == 0);

	# Randomize if this individual shall belong to the population 
        my $rand = rand();
        if ($rand <= ($num_ind_remaining / $num_remaining)) {
                $num_ind_remaining--;

		my @a = split "\t", $l;
		my @g = split " ", $a[5];

		my $IND_ID = $a[0];

		my $phenotype = 0;

		## Add up effects
		foreach my $qtl_id (keys %QTL_LOC) {
			## add genetic effect
			$phenotype += $QTL_EFF{$qtl_id}{$g[$QTL_LOC{$qtl_id}]};
		}

		$ID2GENOPHENO{$IND_ID} = $phenotype;
		$SUM_PHENOTYPES += $phenotype;
		$SQUARED_SUM_PHENOTYPES += $phenotype**2;
	}

	$num_remaining--;

}
close FILE;

print STDERR "Phenotypes of ", keys(%ID2GENOPHENO)+0, " plants recorded.\n";

#########################################################################
## Based on the variance of the phenotypes the environmental
## effect to the phenotype is calculated (why?)

my %ID2ENVPHENO = ();
my %PHENO2ID = ();

# calc mean of phenotypes
my $pheno_mean = $SUM_PHENOTYPES/$IND_COUNT;

# calc variance
my $s2 = (1/($IND_COUNT-1)) * ($SQUARED_SUM_PHENOTYPES - (1/$IND_COUNT) * ($SUM_PHENOTYPES**2) ) ;
my $env_s2 = ($s2 * (1-$HERITABILITY)) / $HERITABILITY;
my $env_s = sqrt($env_s2);

print STDERR "GenoPhenotype mean ", $pheno_mean, " heritability ", $HERITABILITY, " var(geno):", $s2, " var(env):", $env_s2, "\n";

# assess the actual impact
foreach my $ind (keys %ID2GENOPHENO) {
	$ID2ENVPHENO{$ind} = random_normal(1, $pheno_mean, $env_s);
	my $comb_pheno = $ID2ENVPHENO{$ind} + $ID2GENOPHENO{$ind};
	if (defined($PHENO2ID{$comb_pheno})) {
		$PHENO2ID{$comb_pheno} .= "#".$ind;
	}
	else {
		$PHENO2ID{$comb_pheno} = $ind;
	}
}



##############################################################################
# Phenotypic selection of the extreme pools

my $NUM_EXTREME = $IND_COUNT * $EXTREME_PERCENTILE;
print "Number of individuals in each extreme pool: $NUM_EXTREME\n";

# select the individual ids that shall be printed in a second round of file parsing later
my %HIGH_ID = ();
my %LOW_ID = ();

my $count = 0;
LOWO: foreach my $pheno (sort {$a <=> $b} keys %PHENO2ID) { 
	my @ids = split "\n", $PHENO2ID{$pheno};
	for (my $i = 0; $i < @ids; $i++) {
		$LOW_ID{$ids[$i]} = 1;
		$count++;
		if ($count >= $NUM_EXTREME) {
			last LOWO;
		}
	}	
}

$count = 0;
HIGHO: foreach my $pheno (sort {$b <=> $a} keys %PHENO2ID) { 
	my @ids = split "\n", $PHENO2ID{$pheno};
        for (my $i = 0; $i < @ids; $i++) {
	        $HIGH_ID{$ids[$i]} = 1;
		$count++;
        	if ($count >= $NUM_EXTREME) {
                	last HIGHO;
		}
        }
}

#####################################################################################
#Now that the individual ids to be selected are stored parse file again and print respective lines

open OUTH, "> $OUTPUT_FOLDER/Sim.$RUN_ID.high.txt";
open OUTHP, "> $OUTPUT_FOLDER/Sim.$RUN_ID.high.pheno.txt";
open OUTL, "> $OUTPUT_FOLDER/Sim.$RUN_ID.low.txt";
open OUTLP, "> $OUTPUT_FOLDER/Sim.$RUN_ID.low.pheno.txt";

open FILE, $POPULATION_FILE or die "cannot open file $POPULATION_FILE\n";
my $header = <FILE>;
my $cou = 0;
while (my $l = <FILE>) {
	$cou++;
	print $cou, "\n" if $cou%100==0;
	my ($individual) = split " ", $l;
	#my $individual = $a[0];
	if (defined($HIGH_ID{$individual})) {
		print OUTH $l;
		print OUTHP ($ID2GENOPHENO{$individual} + $ID2ENVPHENO{$individual}), "\t", $ID2GENOPHENO{$individual}, "\t", $ID2ENVPHENO{$individual}, "\n";
	}
        if (defined($LOW_ID{$individual})) {
                print OUTL $l;
		print OUTLP ($ID2GENOPHENO{$individual} + $ID2ENVPHENO{$individual}), "\t", $ID2GENOPHENO{$individual}, "\t", $ID2ENVPHENO{$individual}, "\n";
        }
}


open OUTALLPHENO, "> $OUTPUT_FOLDER/Sim.$RUN_ID.all.pheno.txt";
foreach my $key (keys %ID2GENOPHENO) {
	print OUTALLPHENO ($ID2GENOPHENO{$key} + $ID2ENVPHENO{$key}), "\t", $ID2GENOPHENO{$key}, "\t", $ID2ENVPHENO{$key}, "\n";
}
close OUTALLPHENO;


exit(0);

##############################################################################################
# Helper routines

#sub GetCom {
        #my @usage = ("$0 --mapblt file\n\n");

        #die(@usage) if (@ARGV == 0);
        #GetOptions(\%CMD, "mapblt=s");

        #die("Please specify blat file\n") unless defined($CMD{mapblt});

        #$mapblt = $CMD{mapblt};

#}































