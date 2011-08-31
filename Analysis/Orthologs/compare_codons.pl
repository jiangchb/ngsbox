#! /usr/bin/perl

use strict;
use warnings;

use Bio::Perl;
use DBI;

### Connect to A. thaliana and A. lyrata databases
my $dbh_ATH;
my $dbh_Alyr;
&connect_to_ATH();
&connect_to_Alyr();

### General variables;
my $gap_penalty = 2;
my $mm_penalty = 1;

### Get all files from working directory
my @files = glob("*.1.cds.*");


### Set up 4fold degen sites by codon
my %FFD = ();
my @nucs = ("A", "C", "G", "T");

for (my $first = 0; $first <= 3; $first++) {
	my $seq = $nucs[$first];
	for (my $sec = 0; $sec <= 3; $sec++) {
		my $seq2 = $seq . $nucs[$sec];
		for (my $third = 0; $third <= 3; $third++) {	
			my $codon = $seq2 . $nucs[$third];

			my $pep = Bio::Perl::translate_as_string($codon);

			#print STDERR $codon, "\t", $pep, "\n";		
	
			### Check each position for exchanged nucs
			my $codon_string = "";
			for (my $i = 0; $i <= 2; $i++) {
				#print STDERR "   Test ", $i+1, " base:\n";
				my $fourfolddeg = 1;
				for (my $test_nuc = 0; $test_nuc <= 3; $test_nuc++) { 
					my $test_codon = substr($codon, 0, $i) . $nucs[$test_nuc] . substr($codon, $i+1, 3-$i);
					my $test_pep = Bio::Perl::translate_as_string($test_codon);
					#print STDERR "   ", $test_codon, "\t", $test_pep, "\n";
					if ($pep ne $test_pep) {
						$fourfolddeg = 0;
					}
				}
				if ($fourfolddeg == 1) {
					#print STDERR "FOUFOLDDEGENERATE\n";
					$codon_string .= "1";
				}
				else {
					$codon_string .= "0";
				}
			}

			if ($codon_string =~ m/1/) {
				$FFD{$codon} = $codon_string;
				#print STDERR $codon, "\t", $pep, "\t", $codon_string, "\n";		
			}

		}
	}
}


### Find 4fold degen sites by ortholog
CDS: foreach my $file (@files) {

	### Get ATG
	my $agi = substr($file, 0, 9);
	#print STDERR $agi, "\n";
	
	### Set correct isoform
	my $isoform = 1;
	my $best_fit_isoform = 1;
	my $min_edit_distance = 100000;	

	while (-e "$agi.$isoform.cds.fasta") {

		open FILE, "$agi.$isoform.cds.fasta";
		my $ath = "";
		my $lyr = "";
		my $switch = 0;
		while (my $line = <FILE>) {
			chomp($line);
			if (substr($line, 0, 1) eq ">") {
				$switch++;
				if ($line eq ">fake_ortholog") {
					#print STDERR "fake_ortholog\n";
					next CDS;
				}
			}
			else {
				if ($switch == 1) {
					$ath .= $line;
				}
				else {
					$lyr .= $line;
				}
			}
		}
		close FILE;

		### Parse Alignment
		my $edit_distance = 0;
		for (my $i = 0; $i<length($ath); $i++) {
			if (substr($ath, $i, 1) ne substr($lyr, $i, 1)) {
				if (substr($ath, $i, 1) eq "-" or substr($lyr, $i, 1) eq "-") {
					$edit_distance += $gap_penalty;
				}
				else {
					$edit_distance += $mm_penalty;
				}
			}
		}
		

		if ($edit_distance < $min_edit_distance) {
			$min_edit_distance = $edit_distance;
			$best_fit_isoform = $isoform;
		}

		$isoform++;
	}


	### Read in correct alignments		
	my $ath = "";
	my $lyr = "";
	my $switch = 0;
	my $lgi = "";

	open FILE, "$agi.$best_fit_isoform.cds.fasta" or die "Does not exist\n";

	#print STDERR "  $agi.$best_fit_isoform.cds.fasta\n";

	while (my $line = <FILE>) {
		chomp($line);
		if (substr($line, 0, 1) eq ">") {
			if ($switch == 1) {
				$lgi = substr($line, 1, length($line)-1);
			}
			$switch++;
		}
		else {
			if ($switch == 1) {
				$ath .= $line;
			}
			else {
				$lyr .= $line;
			}
		}
	}


	### Parse ALignments
	my $bad_alignment_counter = 0;
	my $ath_codon_pos = 0;
	my $lyr_codon_pos = 0;
	for (my $i = 0; $i<length($ath); $i+=3) {

		my $ath_codon = substr($ath, $i, 3);
		my $lyr_codon = substr($lyr, $i, 3);

		### Check current codon for 4fold degenerate sites
		if (($ath_codon !~ m/-/ || $ath_codon eq "---") and ($lyr_codon !~ m/-/ || $lyr_codon eq "---") and (!($ath_codon =~ m/[^ACGT]/)) and (!($lyr_codon =~ /[^ACGT]/)) ) {

			if ($ath_codon ne "---" and $lyr_codon ne "---") {
				if(	( Bio::Perl::translate_as_string($ath_codon) eq Bio::Perl::translate_as_string($lyr_codon) ) and 
					( defined($FFD{$ath_codon})  and defined($FFD{$lyr_codon}) )
				){
					&get_genome_locus($agi, $best_fit_isoform, $ath_codon_pos+2, $ath_codon, $lgi, $lyr_codon_pos+2, $lyr_codon);
				}	
			}
			
		}
		#else { print STDERR "BADBADBAD:", $agi, " ", $best_fit_isoform, "\n"; exit(1); }

		### Positions
		if ($ath_codon ne "---") {
			$ath_codon_pos += 3;
		}
		if ($lyr_codon ne "---") {
                        $lyr_codon_pos += 3;
                }
	}
}


sub get_genome_locus
{
	my ($agi, $isoform, $ath_cds_pos, $ath_codon, $lgi, $alyr_cds_pos, $alyr_codon) = @_;

	### Arabidopsis thaliana annotation -----------------------------------------------------------------------------
	
	my $ath_chr = -1;
	my $ath_genomic_pos = -1;
	my $ath_ori = "";
	my $ath_snp_change = "";
	my @ath_exons = ();

	### Get A. thaliana exons
	my $q= "SELECT chromosome, begin, end, orientation
		FROM ann_T8_gene
		WHERE tair_id = '$agi' && isoform = $isoform && segment_type = 'CDSexon'
		ORDER BY begin";
	
	my $sth = $dbh_ATH->prepare($q);
	$sth->execute();

	while( my $ref = $sth->fetchrow_hashref() ) {
		$ath_chr = $ref->{chromosome};
		$ath_ori = $ref->{orientation};
		push( @ath_exons, $ref->{begin} );
		push( @ath_exons, $ref->{end} );
	}

	### Calculate genomic position for genes on forward strand
	if( $ath_ori eq "+" ) {
		my $current_cds_pos = 0;
		$ath_genomic_pos = $ath_exons[0];

		for(my $i = 0; $i < @ath_exons; $i +=2) {
			my $exon_len = $ath_exons[$i+1] - $ath_exons[$i] + 1;
			$current_cds_pos += $exon_len;
			if( $i > 0 ) { 
				$ath_genomic_pos += ($ath_exons[$i] - $ath_exons[$i-1] + $exon_len - 1); 
			}
			else {
				$ath_genomic_pos += $exon_len;
			}

			if( $current_cds_pos >= $ath_cds_pos ) {
				$ath_genomic_pos -= ($current_cds_pos - $ath_cds_pos);
				last;
			}
		}
	}

	### Calculate genomic position for genes on reverse strand
	elsif( $ath_ori eq "-" ) {
		my $current_cds_pos = 0;
		$ath_genomic_pos = $ath_exons[@ath_exons-1];

		for(my $i = @ath_exons-1; $i >= 0; $i -=2) {
			my $exon_len = $ath_exons[$i] - $ath_exons[$i-1] + 1;
			$current_cds_pos += $exon_len;
			if( $i < @ath_exons-1 ) {
				$ath_genomic_pos -= ($ath_exons[$i+1] - $ath_exons[$i] + $exon_len - 1);
			}
			else {
				$ath_genomic_pos -= $exon_len;
			}
		
			if( $current_cds_pos >= $ath_cds_pos ) {
				$ath_genomic_pos += ($current_cds_pos - $ath_cds_pos);
				last;
			}
		}
	}

	### Find SNPs in A. thaliana strains at the 4fold degenerate site
	$q = "  SELECT ecotype, reference, basecall FROM poly_snp where chromosome = $ath_chr && position = $ath_genomic_pos && (ecotype = 'Bur-0' || ecotype = 'Tsu-1')"; 
	$sth = $dbh_ATH->prepare($q);
	$sth->execute();
	while( my $ref_ath = $sth->fetchrow_hashref() ) {
		$ath_snp_change .= $ref_ath->{ecotype} .":". $ref_ath->{reference} ."->". $ref_ath->{basecall} .",";
	}
	chop($ath_snp_change);
	if($ath_snp_change eq "") { $ath_snp_change = "NA"; }


	### Arabidopsis lyrata annotation ------------------------------------------------------------------------------------
	
	my $alyr_chr = -1;
	my $alyr_genomic_pos = -1;
	my $alyr_ori = "";
	my $alyr_snp_change = "NA";
	my $alyr_snp_freq = 0;
	my @alyr_exons;

	### Get A. lyrata exons
	$q = "	SELECT chromosome, begin, end, orientation
		FROM Araly1_GeneModels_FilteredModels6
		WHERE id = '$lgi' && isoform = 1 && segment_type = 'CDSexon'
		ORDER BY begin";
	$sth = $dbh_Alyr->prepare($q);
	$sth->execute();

	while( my $ref = $sth->fetchrow_hashref() ) {
		$alyr_chr = $ref->{chromosome};
		$alyr_ori = $ref->{orientation};
		push( @alyr_exons, $ref->{begin} );
		push( @alyr_exons, $ref->{end} );
	}

	### Calculate genomic position for genes on forward strand
	if( $alyr_ori eq "+" ) {
		my $current_cds_pos = 0;
		$alyr_genomic_pos = $alyr_exons[0];

		for(my $i = 0; $i < @alyr_exons; $i +=2) {
			my $exon_len = $alyr_exons[$i+1] - $alyr_exons[$i] + 1;
			$current_cds_pos += $exon_len; 
			if( $i > 0 ) { 
				$alyr_genomic_pos += ($alyr_exons[$i] - $alyr_exons[$i-1] + $exon_len - 1);
			}
			else {  
				$alyr_genomic_pos += $exon_len;
			}

			if( $current_cds_pos >= $alyr_cds_pos ) {
				$alyr_genomic_pos -= ($current_cds_pos - $alyr_cds_pos);
				last;
			}
		}
	}

	### Calculate genomic position for genes on reverse strand
	elsif( $alyr_ori eq "-" ) { 
		my $current_cds_pos = 0;
		$alyr_genomic_pos = $alyr_exons[@alyr_exons-1];

		for(my $i = @alyr_exons-1; $i >= 0; $i -=2) {
			my $exon_len = $alyr_exons[$i] - $alyr_exons[$i-1] + 1;
			$current_cds_pos += $exon_len;
			if( $i < @alyr_exons-1 ) {
				$alyr_genomic_pos -= ($alyr_exons[$i+1] - $alyr_exons[$i] + $exon_len - 1);
			}
			else {
				$alyr_genomic_pos -= $exon_len;
			}

			if( $current_cds_pos >= $alyr_cds_pos ) {
				$alyr_genomic_pos += ($current_cds_pos - $alyr_cds_pos);
				last;
			}
		}
	}

	### Find SNPs in A. lyrata strains at the 4fold degenerate site
	$q = "  SELECT ref_base, major_base, major_freq, minor_base, minor_freq FROM EU_allele_frequency where scaffold = $alyr_chr && position = $alyr_genomic_pos";       
	$sth = $dbh_Alyr->prepare($q);
	$sth->execute();
	my $ref_alyr = $sth->fetchrow_hashref();
	if( exists $ref_alyr->{ref_base} ) {
		if($ref_alyr->{ref_base} ne $ref_alyr->{major_base}) {
			$alyr_snp_change = $ref_alyr->{ref_base} . "->". $ref_alyr->{major_base};
			$alyr_snp_freq = $ref_alyr->{major_freq};
		}
		else {
			$alyr_snp_change = $ref_alyr->{ref_base} . "->". $ref_alyr->{minor_base};
			$alyr_snp_freq = $ref_alyr->{minor_freq};
		}
	}


	### -------------------------------------------------------------------------------------------------------------
	### Print:
	### AGI, iso, ori, cds pos, chr, genomic pos, snp count, Ath codon seq, 
	### LGI, iso, ori, cds pos, chr, genomic pos, snp count, Alyr codon seq
	print 	"$agi\t$isoform\t$ath_ori\t$ath_cds_pos\t$ath_chr\t$ath_genomic_pos\t$ath_snp_change\t$ath_codon\t" .
		"$lgi\t1\t$alyr_ori\t$alyr_cds_pos\t$alyr_chr\t$alyr_genomic_pos\t$alyr_snp_change\t$alyr_snp_freq\t$alyr_codon\n";
}


############################
### Connects to ATH database
sub connect_to_ATH 
{
	my $databaseName = "solexa";
	my $driver = "mysql";
	my $host = "orb.eb.local";
	my $username = "solexa";
	my $password = "s0lexa";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh_ATH = DBI->connect($dsn, $username, $password ) || db_connect_error();
}

#############################
### Connects to Alyr database
sub connect_to_Alyr 
{
	my $databaseName = "Alyr_Pool";
	my $driver = "mysql";
	my $host = "orb.eb.local";
	my $username = "solexa";
	my $password = "s0lexa";
	my $dsn = "DBI:$driver:database=$databaseName;host=$host";
	my $drh = DBI->install_driver("mysql");
	$dbh_Alyr = DBI->connect($dsn, $username, $password ) || db_connect_error();
}







