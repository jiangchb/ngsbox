#! /usr/bin/perl
use strict;
use lib "$ENV{PGSP}/Assembly/Realign/";
use needle_alignment;

open CONTIG, $ARGV[0];
open SANGERSHORT, $ARGV[1] or die "C";

my $count_seq = 0;
my $count_perfects = 0;
my $count_nearly_perfect = 0;

my %SEQUENCE = ();

my $id = "";
my $seq = "";

while (my $line = <CONTIG>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		if ($id ne "") {
			$SEQUENCE{$id} = $seq;
		}
		$id = substr($line, 1, length($line)-1);
		$seq = "";
	}
	else {
		$seq .= $line;
	}
	
}
if ($id ne "") {
	$SEQUENCE{$id} = $seq;
}


$id = "";
$seq = "";


while (my $line = <SANGERSHORT>) {
	chomp($line);
        if (substr($line, 0, 1) eq ">") {
		if ($id ne "") {
			align($id, $SEQUENCE{$id}, $seq);	
		}
                $id = substr($line, 1, length($line)-1);
                $seq = "";
	}
	else {
                $seq .= $line;
        }
}
if ($id ne "") {
	align($id, $SEQUENCE{$id}, $seq);
}


print "Sequences: ", $count_seq, "\n";
print "Identical: ", $count_perfects, "\n";
print "Nearly identical: ", $count_nearly_perfect, "\n";


sub align {
	my ($id, $seq1, $seq2) = @_;


	my $aligner = new needle_alignment();
	$aligner->needleman_wunsch($seq1, $seq2);
	$aligner->parse_fasta_alignment();

	#my $align = needle_alignment::needleman_wunsch($seq1, $seq2);
	#my ($snp_ref, $del_ref, $ins_ref) = needle_alignment::parse_alignment($align, 0);

	if (not ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0)) {
		print "\n\n";
		print ">$id\n";
		print $aligner->{align_seq1}, "\n";
        	print $aligner->{align_seq2}, "\n";
		for (my $i = 0; $i<length($aligner->{align_seq1}); $i++) {
			if (uc(substr($aligner->{align_seq2}, $i, 1)) ne uc(substr($aligner->{align_seq1}, $i, 1))) {
				print "X";
			}
			else {
				print " ";
			}
		}
		print "\n";

		print "NumSNPS:", $aligner->{num_subs}, "\tNumDels:", $aligner->{num_del}, "\tNumIns:", $aligner->{num_ins}, "\n";
	}


	$count_seq++;		
	if ($aligner->{num_subs} == 0 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0) {
		$count_perfects++;
	}
	elsif ($aligner->{num_subs} < 3 and $aligner->{num_del} == 0 and  $aligner->{num_ins} == 0) {
		$count_nearly_perfect++;
	}

}


