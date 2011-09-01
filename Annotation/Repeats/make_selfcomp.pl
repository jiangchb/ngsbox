#!/usr/bin/perl
####################################################################################
#Author 	Korbinian Schneeberger 
#Date 		05/15/07
#Version	0.1
#Input		Fasta mit Genome reads	
#Function	Performs mapping and outputs h0, h1, h2, h3 repeat mappings
####################################################################################

use strict;
use warnings;
use Getopt::Long;
use Cwd;

my %CMD = ();
my $len;
my $idx;
my $mapfile;
my $vmatchfile;
my $map_flag = 0;
my $seg;
my $prb_string;
my $max_mm;

GetCom();

##########
# Perform Mapping
##########

if ($map_flag ==  1) {	
	if ($max_mm > 0) {
		system("vmatch -q $mapfile -d -p -l ".$len." -h $max_mm -s abbrev -showdesc 30 -seedlength 9 -noscore -noidentity $idx | grep '.' > $vmatchfile");
	}
	else {
		system("vmatch -q $mapfile -d -p -l ".$len." -s abbrev -showdesc 30 -noscore -noidentity $idx | grep '.' > $vmatchfile");
	}

	clean_vmatch($vmatchfile);
	system("sort --buffer-size=10% -k6 $vmatchfile.parsed > tmp");
	system("mv tmp $vmatchfile.parsed");
	$vmatchfile .= ".parsed"; 
}

open FILE, "$vmatchfile" or die "Cannot open file $vmatchfile\n";

open MAPPING_0, "> map.mm0.list.unsorted" or die "cannot open file\n";

if ($seg == 0) {
	if ($max_mm > 0) {
		open MAPPING_1, "> map.mm1.list.unsorted" or die "cannot open file\n";
	}
	if ($max_mm > 1) {
		open MAPPING_2, "> map.mm2.list.unsorted" or die "cannot open file\n";
	} 
	if ($max_mm > 2) {
		open MAPPING_3, "> map.mm3.list.unsorted" or die "cannot open file\n";
	}
} else {
	$prb_string = make_prb_string($len); 
}

my $curr_frag = "-1";

my $MIN_EVALUE_0;
my $MIN_EVALUE_1;
my $MIN_EVALUE_2;
my $MIN_EVALUE_3;

my @FRAG_ID_0;
my @ALIGN_0;
my @CHR_0;
my @CHRPOS_0;
my @ORI_0;
my @LEN_0;
my @OFFSET_0;

my @FRAG_ID_1;
my @ALIGN_1;
my @CHR_1;
my @CHRPOS_1;
my @ORI_1;
my @LEN_1;
my @OFFSET_1;

my @FRAG_ID_2;
my @ALIGN_2;
my @CHR_2;
my @CHRPOS_2;
my @ORI_2;
my @LEN_2;
my @OFFSET_2;

my @FRAG_ID_3;
my @ALIGN_3;
my @CHR_3;
my @CHRPOS_3;
my @ORI_3;
my @LEN_3;
my @OFFSET_3;

while (<FILE>) {
	my @val = split " ", $_;
	
	my $length = $val[0];
	my $chr = $val[1];
	my $chrpos = $val[2];
	my $orientation = $val[3];
	my $frag = $val[5];
	my $offset = $val[6];
	my $mm = $val[7];
	my $evalue = $val[8];
	my ($frag_chr, $frag_pos) = split "_", $frag;
	$chrpos = $chrpos + 1;
	my $align = $val[9];

#print STDERR "0:",$len, "\n";
#"\t1:",$length, "\t2:", $chr, "\t3:", $chrpos, "\t4:", $orientation, "\t5:", $frag, "\t6:", $mm, "\t7:", $evalue, "\t8:", $align, "\n";

	if ((not ($seg == 0 and $chr == $frag_chr and $chrpos == $frag_pos)) and $length == $len) {

		if ($curr_frag ne $frag) {
			if ($curr_frag ne "-1") {
				print_mapping($max_mm);
			}
			reset_evalues();
			reset_values(\@FRAG_ID_0, \@ALIGN_0, \@CHR_0, \@CHRPOS_0, \@ORI_0, \@LEN_0, \@OFFSET_0);
			reset_values(\@FRAG_ID_1, \@ALIGN_1, \@CHR_1, \@CHRPOS_1, \@ORI_1, \@LEN_1, \@OFFSET_1);
			reset_values(\@FRAG_ID_2, \@ALIGN_2, \@CHR_2, \@CHRPOS_2, \@ORI_2, \@LEN_2, \@OFFSET_2);
			reset_values(\@FRAG_ID_3, \@ALIGN_3, \@CHR_3, \@CHRPOS_3, \@ORI_3, \@LEN_3, \@OFFSET_3);
			$curr_frag = $frag;
		} 

		if (abs($mm) == 0) {
			if ($evalue < $MIN_EVALUE_0) {
				reset_values(\@FRAG_ID_0, \@ALIGN_0, \@CHR_0, \@CHRPOS_0, \@ORI_0, \@LEN_0, \@OFFSET_0);
			} 
			if ($evalue <= $MIN_EVALUE_0) {
				set_values($frag, $align, $chr, $chrpos, $orientation, $length, $offset, \@FRAG_ID_0, \@ALIGN_0, \@CHR_0, \@CHRPOS_0, \@ORI_0, \@LEN_0, \@OFFSET_0);
				$MIN_EVALUE_0 = $evalue;
			}
		} 
		if ($seg == 0) {
			if (abs($mm) == 1) {
				if ($evalue < $MIN_EVALUE_1) {
                		        reset_values(\@FRAG_ID_1, \@ALIGN_1, \@CHR_1, \@CHRPOS_1, \@ORI_1, \@LEN_1, \@OFFSET_1);
		                } 
				if ($evalue <= $MIN_EVALUE_1) {
                		        set_values($frag, $align, $chr, $chrpos, $orientation, $length, $offset, \@FRAG_ID_1, \@ALIGN_1, \@CHR_1, \@CHRPOS_1, \@ORI_1, \@LEN_1, \@OFFSET_1);
					$MIN_EVALUE_1 = $evalue;
		                }
			} 
			if (abs($mm) == 2) {
				if ($evalue < $MIN_EVALUE_2) {
                	        	reset_values(\@FRAG_ID_2, \@ALIGN_2, \@CHR_2, \@CHRPOS_2, \@ORI_2, \@LEN_2, \@OFFSET_2);
		                } 
				if ($evalue <= $MIN_EVALUE_2) {
                		        set_values($frag, $align, $chr, $chrpos, $orientation, $length, $offset, \@FRAG_ID_2, \@ALIGN_2, \@CHR_2, \@CHRPOS_2, \@ORI_2, \@LEN_2, \@OFFSET_2);
					$MIN_EVALUE_2 = $evalue;
		                }
        		} 
			if (abs($mm) == 3) {
				if ($evalue < $MIN_EVALUE_3) {
                	        	reset_values(\@FRAG_ID_3, \@ALIGN_3, \@CHR_3, \@CHRPOS_3, \@ORI_3, \@LEN_3, \@OFFSET_3);
		                } 
				if ($evalue <= $MIN_EVALUE_3) {
                		        set_values($frag, $align, $chr, $chrpos, $orientation, $length, $offset, \@FRAG_ID_3, \@ALIGN_3, \@CHR_3, \@CHRPOS_3, \@ORI_3, \@LEN_3, \@OFFSET_3);
					$MIN_EVALUE_3 = $evalue;
		                }
        		}
		}
	}

}
print_mapping($max_mm);

if ($seg == 0) {
	system("sort --buffer-size=10% -k1n -k2n map.mm0.list.unsorted > map.mm0.list.sorted");
	system("rm map.mm0.list.unsorted");
}
else {
	system("sort --buffer-size=10% -k1n -k2n map.mm0.list.unsorted > map.mm0.list.segmentation.sorted");
        system("rm map.mm0.list.unsorted");
}

if ($seg == 0) {
	if ($max_mm > 0) {
		system("sort --buffer-size=10% -k1n -k2n map.mm1.list.unsorted > map.mm1.list.sorted");
		system("rm map.mm1.list.unsorted");
	}
	if ($max_mm > 1) {
		system("sort --buffer-size=10% -k1n -k2n map.mm2.list.unsorted > map.mm2.list.sorted");
		system("rm map.mm2.list.unsorted");
	}
	if ($max_mm > 2) {
		system("sort --buffer-size=10% -k1n -k2n map.mm3.list.unsorted > map.mm3.list.sorted");
		system("rm map.mm3.list.unsorted");
	}

	if ($max_mm > 0) {
		system("sort --buffer-size=10% -m -k1n -k2n map.mm0.list.sorted map.mm1.list.sorted > map.mm0_1.list.sorted");
	}
	if ($max_mm > 1) {
		system("sort --buffer-size=10% -m -k1n -k2n map.mm0.list.sorted map.mm1.list.sorted map.mm2.list.sorted > map.mm0_1_2.list.sorted");
	}
	if ($max_mm > 2) {
		system("sort --buffer-size=10% -m -k1n -k2n map.mm0.list.sorted map.mm1.list.sorted map.mm2.list.sorted map.mm3.list.sorted > map.mm0_1_2_3.list.sorted");
	}

}


exit(0);

sub print_mapping {
	my ($max_mm) = @_;

	if ($max_mm >= 0) {
		for (my $i = 0; $i<@FRAG_ID_0; $i++) {
			print MAPPING_0 $CHR_0[$i], "\t", $CHRPOS_0[$i], "\t", $ALIGN_0[$i], "\t", $FRAG_ID_0[$i], "\t", $ORI_0[$i], "\t", "0", "\t", @FRAG_ID_0 + 0, "\t", $LEN_0[$i], "\t", $OFFSET_0[$i];
			if ($seg == 1) {
				print MAPPING_0 "\t", $prb_string, "\n"; 
			} else {
				print MAPPING_0 "\n";
			}
		}	
	}

	if ($max_mm > 0) {
		for (my $i = 0; $i<@FRAG_ID_1; $i++) {
			print MAPPING_1 $CHR_1[$i], "\t", $CHRPOS_1[$i], "\t", $ALIGN_1[$i], "\t", $FRAG_ID_1[$i], "\t", $ORI_1[$i], "\t", "1", "\t", @FRAG_ID_1 + 1, "\t", $LEN_1[$i], "\t", $OFFSET_1[$i], "\n";
		}
	}

	if ($max_mm > 1) {
		for (my $i = 0; $i<@FRAG_ID_2; $i++) {
			print MAPPING_2 $CHR_2[$i], "\t", $CHRPOS_2[$i], "\t", $ALIGN_2[$i], "\t", $FRAG_ID_2[$i], "\t", $ORI_2[$i], "\t", "2", "\t", @FRAG_ID_2 + 2, "\t", $LEN_2[$i], "\t", $OFFSET_2[$i], "\n";
		}
	}

	if ($max_mm > 2) {
		for (my $i = 0; $i<@FRAG_ID_3; $i++) {
			print MAPPING_3 $CHR_3[$i], "\t", $CHRPOS_3[$i], "\t", $ALIGN_3[$i], "\t", $FRAG_ID_3[$i], "\t", $ORI_3[$i], "\t", "3", "\t", @FRAG_ID_3 + 3, "\t", $LEN_3[$i], "\t", $OFFSET_3[$i], "\n";
		}
	}
}

sub set_values {
	my ($val1, $val2, $val3, $val4, $val5, $val6, $val7, $ref1, $ref2, $ref3, $ref4, $ref5, $ref6, $ref7) = @_;

	push @$ref1, $val1;
	push @$ref2, $val2;
	push @$ref3, $val3;
	push @$ref4, $val4;
	push @$ref5, $val5;
	push @$ref6, $val6;
	push @$ref7, $val7;

	return(0);
}

sub reset_evalues {

	$MIN_EVALUE_0 = 10;
	$MIN_EVALUE_1 = 10;
	$MIN_EVALUE_2 = 10;
	$MIN_EVALUE_3 = 10;

	return(0);
}

sub reset_values {

	foreach my $ref (@_) {
		@$ref = ();
	}

	return(0);
}

sub clean_vmatch {
	my ($file) = @_;
	
	open FILE, $file or die "Cannot open file\n";
	open OUT, "> $file.parsed";
	
	while (<FILE>) {
		chomp;
		if (substr($_, 0, 1) ne "#") {
			print OUT $_, "\t";
			my $a = <FILE>;
			print OUT $a;
		}
	}

	close FILE or die "Cannot close file\n";
	close OUT or die "Cannot close file\n";

	return(0);
}

sub make_prb_string {
	my ($len) = @_;
	my $s = "S";
	for (my $i = 0; $i<$len-1; $i++) {
		$s .= "S";
	}
	return $s."\t".$s."\t".$s;
}

sub GetCom {

  my @usage = ("\nUsage: $0 

required:
--file\t\tVmatch output or fragments for mapping
--len\t\tHit length (usually = readlength) 
--seg\t'0' or '1'. Used for segmentation. Take self hits into account. Drops MM hits.
--max_mm\tMaximum mismatches in the self mapping (0, 1, 2 or 3)

optional:
--map\t\tperform mapping, if \'1\' also specify:
--index\t\tVmatch index file

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "file=s", "index=s", "len=s", "map=s", "seg=s", "max_mm=s");

	die("Please specify an input file\n") unless defined($CMD{file});
	die("Please specify if the run is a segmentation run or a rep scale calc run.\n") unless defined($CMD{seg});
	die("Please specify maximal mismatches\n") unless defined($CMD{max_mm});

	if (not defined($CMD{map}) or $CMD{map} != 1) {
		die("Please specify a fragment length\n") unless defined($CMD{len});
		$map_flag = 0;
		$vmatchfile = $CMD{file};
		$len = $CMD{len};
	} else {
		die("Please specify a fragment length\n") unless defined($CMD{len});
		die("Please specify a vmatch index\n") unless defined($CMD{index});
		$mapfile = $CMD{file};
		$vmatchfile = $mapfile.".out";
		$map_flag = 1;
		$len = $CMD{len};
		$idx = $CMD{index};
	}

	$seg = $CMD{seg};
	$max_mm = $CMD{max_mm};

	if ($max_mm != 0 and $max_mm != 1 and $max_mm != 2 and $max_mm != 3) {
		die ("Max mismatches must equal 0, 1, 2 or 3\n");
	}

	return(0);
}



exit(0);
