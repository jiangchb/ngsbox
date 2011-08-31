#! /usr/bin/perl

my @ACC = ("Col-0","Agu-1","Bak-2","Bak-7","Cdm-0","Del-10","Dog-4","Don-0","Ey15-2","Fei-0","HKT2.4","ICE1","ICE102","ICE104","ICE106","ICE107","ICE111","ICE112","ICE119","ICE120","ICE127","ICE130","ICE134","ICE138","ICE150","ICE152","ICE153","ICE163","ICE169","ICE173","ICE181","ICE21","ICE212","ICE213","ICE216","ICE226","ICE228","ICE29","ICE33","ICE36","ICE49","ICE50","ICE60","ICE61","ICE63","ICE7","ICE70","ICE71","ICE72","ICE73","ICE75","ICE79","ICE91","ICE92","ICE93","ICE97","ICE98","Istisu-1","Kastel-1","Koch-1","Lag2.2","Leo-1","Lerik1-3","Mer-6","Nemrut-1","Nie1-2","Ped-0","Pra-6","Qui-0","Rue3-1-31","Sha","Star-8","TueSB30-3","Tuescha-9","TueV13","TueWa1-2","Vash-1","Vie-0","WalhaesB4","Xan-1","Yeg-1");

my $usage = "$0 positions_file matrix\n";
my $pos_file = shift or die $usage;
my $matrix_file = shift or die $usage;


open POS, $pos_file or die "Cannot open $pos_file\n";

my %REG = ();
my %IDREG = ();

my $c = 1;
while (<POS>) {
	my @a = split " ";
	for (my $i = $a[1]; $i <= $a[2]; $i++) {
		$REG{$a[0]."\t".$i} = $c;	
	}
	$IDREG{$c} = $a[0]."\t".$a[1]."\t".$a[2];
	$c++;
}
close POS;

open MATRIX, $matrix_file or die "cannot open file\n";

my %SEQ = ();


MAT: while (my $l = <MATRIX>) {
	my @a = split " ", $l;
	if ($a[1]%100000==0) { 
		print STDERR $a[0], "\t", $a[1], "\n"; 
	}
	
	if (defined($REG{$a[0]."\t".$a[1]})) {
		for (my $i = 2; $i <= 82; $i++) {
			my $char = $a[$i];
			if ($char ne "A" and $char ne "C" and $char ne "G" and $char ne "T") {
				$char = "N";
			}
			$SEQ{$ACC[$i-2]} .= $char;
		}
		delete $REG{$a[0]."\t".$a[1]};
	}

	# exit loop if no region left
	if (keys(%REG) == 0) {
		last MAT;
	}

}

for (my $i = 0; $i < @ACC; $i++) {
	print ">".$ACC[$i]."\n";
	print rev_comp($SEQ{$ACC[$i]}), "\n";
}



sub rev_comp {
        my ($seq) = @_;
        my $new = reverse($seq);
        $new =~ tr/acgtACGT/tgcaTGCA/;

        return $new;
}




