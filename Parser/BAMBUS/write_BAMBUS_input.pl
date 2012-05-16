#! /usr/bin/perl

use strict;

my $usage = "$0 prefix assembly.fa libname lib_min lib_max out1.blat out2.blat read_pair1.fa read_pair2.fa format\n";
if ((@ARGV+0) > 0) {
	write_log();
}

my $prefix = shift or die $usage;
my $assm_file = shift or die $usage;

## foreach lib
my $lib = shift or die $usage;
my $lib_min = shift or die $usage;
my $lib_max = shift or die $usage;
my $align1_file = shift or die $usage;
my $align2_file = shift or die $usage;
my $reads1_file = shift or die $usage;
my $reads2_file = shift or die $usage;

my $format = shift or die $usage;

# sequences
my %READSEQ = ();
my %CONTIGSEQ = ();

# alignments
my %RCOUNT = ();
my %C2R = (); # Contig to read alignments
my %R2C = (); # Read to contig



## foreach lib
read_fasta(\%READSEQ, $lib, "fwd", $reads1_file);
read_fasta(\%READSEQ, $lib, "rev", $reads2_file);
read_fasta(\%CONTIGSEQ, "assembly", "", $assm_file);

read_alignments($format, \%C2R, \%R2C, \%RCOUNT, $lib, "fwd", $align1_file);
read_alignments($format, \%C2R, \%R2C, \%RCOUNT, $lib, "rev",  $align2_file);

write_conf_file();
write_mates_file();
write_contig_file();


exit(0);



sub write_contig_file {
	open OUT, ">$prefix.contig" or die "cannot write contig file\n";

	my $out = "";
	my $flag = 0;

	my $count_contig = 0;
	my $count_contig_hit  = 0; # flag count
	my $count_contig_hits = 0; # flag count

	my $count_alignments = 0;
	my $count_alignments_mapped_pair = 0;
	my $count_alignments_mapped_pair_uniq = 0;
	my $count_alignments_mapped_pair_uniq_diff = 0;

	my $count_contigpairs = 0;

	my %CTG_PAIRS = ();

	foreach my $con (keys %CONTIGSEQ) {
		# print contig
		$out = "";
		$out .= "##$con 1 ".length($CONTIGSEQ{$con})." bases.\n";
		$out .= fixed_length($CONTIGSEQ{$con}, 80);
		$flag = 0;

		$count_contig++;

		# print reads on contig
		foreach my $lib (keys %{$C2R{$con}}) {
			foreach my $read (keys %{$C2R{$con}{$lib}}) {
				$count_alignments++;
				if (defined($RCOUNT{$lib}{$read}{"fwd"}) and defined($RCOUNT{$lib}{$read}{"rev"})) {
					$count_alignments_mapped_pair++;
					if ($RCOUNT{$lib}{$read}{"fwd"} == 1 and $RCOUNT{$lib}{$read}{"rev"} == 1) {
						$count_alignments_mapped_pair_uniq++;
				     		if ($R2C{$lib}{$read}{"fwd"} ne $R2C{$lib}{$read}{"rev"}) {
							$count_alignments_mapped_pair_uniq_diff++;
							foreach my $ori (keys %{$C2R{$con}{$lib}{$read}}) {
								$out .= $C2R{$con}{$lib}{$read}{$ori}."\n";
								$out .= fixed_length($READSEQ{$lib}{$read}{$ori}, 80);
							}
							$flag=1;
							$count_contig_hits++;

							## count bridges and bridge strength

							my $pid = ""; 
							if ($R2C{$lib}{$read}{"fwd"} < $R2C{$lib}{$read}{"rev"}) {
								$pid = $R2C{$lib}{$read}{"fwd"}."#".$R2C{$lib}{$read}{"rev"};
							}
							else {
								$pid = $R2C{$lib}{$read}{"rev"}."#".$R2C{$lib}{$read}{"fdw"};
							}

							if (not defined($CTG_PAIRS{$pid})) {
								$count_contigpairs++;
							}
	
							$CTG_PAIRS{$pid}++;
						}
					}
				}
			}
		}

		if ($flag == 1) {
			$count_contig_hit++;
		}
		
		print OUT $out;

	}

	
	print "\n";
	print "# Contigs: $count_contig\n";
	print "# Contigs with at least one useful mate: $count_contig_hit\n";
	my $avg_hits_contig = $count_contig_hit == 0 ? 0 : ($count_contig_hits/$count_contig_hit);
	print "# Avg mates per contig: $avg_hits_contig\n";
	print "\n";
	print "Contig pairs: $count_contigpairs\n";
	my $avg_bridges_contigpair = ($count_alignments_mapped_pair_uniq_diff/$count_contigpairs);
	print "Avg bridges per contig pairs: $avg_bridges_contigpair\n";
	print "\n";
	print "# Alignments: $count_alignments\n";
	print "# Alignments mapped: $count_alignments_mapped_pair\n";
	print "# Alignments mapped uniq: $count_alignments_mapped_pair_uniq\n";
	print "# Alignments mapped uniq (diff contigs): $count_alignments_mapped_pair_uniq_diff\n";

	close OUT;
	
	open OUT, ">".$prefix.".contig_pair_coverage" or die "cannot open contig_pair_coverage file\n";
	foreach my $cp (keys %CTG_PAIRS) {
		print OUT $CTG_PAIRS{$cp}, "\n";
	}
	close OUT;

}

sub write_mates_file {
	open OUT, ">$prefix.mates" or die "cannot write mates file\n";
	
	print OUT "library\t$lib\t$lib_min\t$lib_max\n";

	foreach my $lib (keys %READSEQ) {
		foreach my $read (keys %{$READSEQ{$lib}}) {
			print OUT $read."_fwd\t".$read."_rev\t$lib\n";
		}
	}

	close OUT;
}

sub write_conf_file {
	open OUT, ">$prefix.conf" or die "cannot write conf file\n";

print OUT "
# Priorities

priority lib_$lib 1
priority ALL 1

# Redundancies
# redundancy lib_some 1

# Global redundancy
redundancy 2

# min group size
mingroupsize 0

";

	close OUT;
}

sub read_alignments { 
	my ($format, $hash_ref, $r2c_ref, $count_ref, $lib, $readori, $file) = @_;

	open FILE, $file or die "cannot open file $file\n";

	my $count_c = 0;

	if ($format eq "blat") {
		while (my $l = <FILE>) {
			## skip header
			if (substr($l, 0, 8) eq "psLayout") {
				for (my $i = 0; $i < 5; $i++) {
					$l = <FILE>;
				}
			}
			## parse 
			my @a = split " ", $l;

			my $rid = $a[9];
			my $rstart = $a[11];
	                my $rend = $a[12];
			my $rlen = $a[10];

			my $cid = $a[13];
			my $cstart = $a[15];
			my $cend = $a[16];
			
			my $match = $a[0];
			my $orientation = $a[8]; 
			
			if ($match / $rlen > 0.5) {
				## Store in .contig file format
				my $oristr = $orientation eq "+" ? "[]" : "[RC]";
				my $rhitstr = $orientation eq "+" ? "{$rstart $rend}" : "{$rend $rstart}";
				my $chitstr = "<$cstart $cend>";
	
				if (not defined(${$hash_ref}{$cid})) {
					$count_c++;
				}

				# Bug? Multiple hits on one contig?
				${$hash_ref}{$cid}{$lib}{$rid}{$readori} = "#$rid"."_"."$readori($cstart) $oristr $rlen bases. $rhitstr $chitstr";
				${$r2c_ref}{$lib}{$rid}{$readori} = $cid; 
				${$count_ref}{$lib}{$rid}{$readori}++;
			}
		} 
	}
	elsif ($format eq "map.list") {

		while (my $l = <FILE>) {
			my @a = split " ", $l;
	
	                my $rid = $a[3];
        	        my $rstart = 1;
                	my $rend = $a[7];
	                my $rlen = $a[7];

	                my $cid = $a[0];
        	        my $cstart = $a[1];
                	my $cend = $a[1] + $a[7];
                        
	                my $match = $a[7] - $a[5];
        	        my $orientation = $a[4] eq "D" ? "+" : "-";

			## Store in .contig file format
        	        my $oristr = $orientation eq "+" ? "[]" : "[RC]";
                	my $rhitstr = $orientation eq "+" ? "{$rstart $rend}" : "{$rend $rstart}";
	                my $chitstr = "<$cstart $cend>";

	
	                if (not defined(${$hash_ref}{$cid})) {
				$count_c++;
                	}

	                # Bug? Multiple hits on one contig?
        	        ${$hash_ref}{$cid}{$lib}{$rid}{$readori} = "#$rid"."_"."$readori($cstart) $oristr $rlen bases. $rhitstr $chitstr";
                	${$r2c_ref}{$lib}{$rid}{$readori} = $cid;
	                ${$count_ref}{$lib}{$rid}{$readori}++;
		}
	}
	else {
		die("alignment format needs to be blat or map.list formated files.\n");
	}

	print "Found $count_c novel contigs\n";

}

sub read_fasta {
	my ($hash_ref, $lib, $id_ext, $file) = @_;

	open FILE, $file or die "cannot open file $file\n";

	my $seq = "";
	my $id = "";
	my $count_e = 0;

	while (<FILE>) {
		chomp($_);
		if (substr($_, 0, 1) eq ">") {
			my @a = split " ", $_;
			if ($seq ne "") {
				$count_e++;
				if ($lib eq "assembly") {
					${$hash_ref}{$id} = $seq;
				}
				else {
					${$hash_ref}{$lib}{$id}{$id_ext} = $seq;
				}
			}
			$id = substr($a[0], 1, length($a[0])-1);
			$seq = "";
		}	
		else {
			$seq .= $_;
		}
	}

	if ($seq ne "") {
		$count_e++;
		if ($lib eq "assembly") {
			${$hash_ref}{$id} = $seq;
                }
                else {
                	${$hash_ref}{$lib}{$id}{$id_ext} = $seq;
                }
	}
	
	print "Read fasta, found $count_e\n";	

}

sub fixed_length {
	my ($str, $len) = @_;
	my $ret = "";
	for (my $i = 1; $i <= length($str); $i++) {
		$ret .= substr($str, $i, 1);
		if ($i % $len == 0) {
			$ret .= "\n";
		}
	}
	if (length($str) % $len != 0) {
		$ret .= "\n";
	}

	return $ret;
}


sub write_log {
	my $logf = $ARGV[0].".log";
	open OUT, ">>$logf" or die "cannot open log file\n";
	print OUT "perl $0 ";
	print OUT join(" ", @ARGV), "\n";
	close OUT;
}


