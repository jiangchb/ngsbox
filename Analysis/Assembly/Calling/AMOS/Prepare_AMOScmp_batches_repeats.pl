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
#  Module: Analysis::Assembly::Calling::AMOS::Prepare_AMOScmp_batches_repeats.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 ReferenceFasta   RepeatFile  LeftOverContigs  OutputFolder(will be created)   AssemblyFolder   [AssemblyFolder [AssemblyFolder [...]]]\n\n";

my $reffasta = shift or die $usage;
my $repeatfile = shift or die $usage;
my $leftoverfile = shift or die $usage;
my $outfolder = shift or die $usage;

my @assemblyfolders = @ARGV;

if (!(-e $outfolder)) {
	mkdir($outfolder);
	mkdir($outfolder."/AMOScmp_batches/");
}
else {
	die "$outfolder already exists\n";
}

######################################################################
open FFILE, $reffasta or die "Cannot open file\n";
my %REFSEQ = ();
my $refid = "";
my $refseq = "";
while (my $line = <FFILE>) {
	chomp($line);
	if (substr($line, 0, 1) eq ">") {
		if ($refid ne "") {
			$REFSEQ{$refid} = $refseq;
		}
		$refseq = "";
		my @a = split " ", substr($line, 1, length($line)-1);
		$refid = $a[0];
	}
	else {
		$refseq .= $line;
	}
}
close FFILE;

######################################################################

open RFILE, $repeatfile or die "Cannot open file\n";
my @REPEAT_CHR = ();
my @REPEAT_BEGIN = ();
my @REPEAT_END = ();
while (<RFILE>) {
	my @a = split " ";
	push @REPEAT_CHR, $a[0];
	push @REPEAT_BEGIN, $a[1];
	push @REPEAT_END, $a[2];
}
close RFILE;

######################################################################

my $diff = 0; my $co = 0;
for (my $asi = 0; $asi < @assemblyfolders-1; $asi++) {
	for (my $asj = $asi+1; $asj < @assemblyfolders; $asj++) {
		$co++;
		my $file1 = $assemblyfolders[$asi]."/superblocks.txt";
		my $file2 = $assemblyfolders[$asj]."/superblocks.txt";
		my $comp = `diff $file1 $file2`;
		if ($comp ne "") {
			print STDERR "[ERROR] Superblock files do not match:\n";
			print STDERR "    ", $file1, "\n    ", $file2, "\n";
			$diff = 1;
		}
	}
}
if ($diff == 0) {
	print STDERR "[INFO] Superblock files match (", $co, " comparisons)\n";
}

######################################################################


my $superblocks = $assemblyfolders[0]."/superblocks.txt";

open FILE, $superblocks or die "Cannot open file: ".$superblocks."\n";

my $last_target_end = -1;

my $batchnum = 1;

my $batch_start = -1;
my $batch_chr = -1;
my @batch_targets = ();

my $target_chr;
my $target_id;
my $target_start;
my $target_end;

my $repeat = 0;


while (<FILE>) {
	my @a = split " ";
	$target_id = $a[0];
	$target_chr = $a[1];
	$target_start = $a[2];
	$target_end = $a[3];
	
	# Set the start of the very first batch
	if ($batch_start == -1) {
		$batch_start = $target_start;
	}

	# The new target does not belong to the current batch: batch gets finalized:
	if ($batch_chr != -1) {
		if ($batch_chr != $target_chr or ($target_end > $REPEAT_BEGIN[$repeat] and $target_chr == $REPEAT_CHR[$repeat])) {
			create_batch($batchnum, \@batch_targets, $batch_chr, $batch_start, $last_target_end);
                	@batch_targets = ();
	                $batch_start = $target_start;
			$repeat++;
			$batchnum++;
		}
		while (($REPEAT_CHR[$repeat] < $target_chr or ($REPEAT_CHR[$repeat] == $target_chr and $REPEAT_BEGIN[$repeat] < $target_end)) and $repeat < @REPEAT_CHR) {
			$repeat++;
		}
	}

	$batch_chr = $target_chr;
	push @batch_targets, $target_id;
	$last_target_end = $target_end;
}

create_batch($batchnum, \@batch_targets, $batch_chr, $batch_start, $last_target_end);



sub create_batch {

	my ($batchnum, $batch_targets_ref, $chr, $begin, $end) = @_;

	my $folder = $outfolder."/AMOScmp_batches/AMOS_batch_".$batchnum;
	if (!(-e $folder)) {
		mkdir($folder);
	}

	open OUT, ">>".$outfolder."/AMOScmp_batches/batch_desc.txt";
	print OUT $chr, "\t", $begin, "\t", $end, "\n";
	close OUT;

	### Collect all contigs belonging to this batch
	my $cmd = "cat ";
	foreach my $assemblyfolder (@assemblyfolders) {
		for (my $i = 0; $i < @$batch_targets_ref; $i++) {

			my $c = 0;

			my $file1 = "$assemblyfolder/BuildingSite/Scaffolds/scaffold_".${$batch_targets_ref}[$i]."/contigs.fa";
			if (-e $file1) {
				$cmd .= " ".$file1;
				$c++;
				if (-s $file1 == 0) {
					print STDERR "[ERROR] ".$file1." is empty\n";
				}
			}

			my $file2 = "$assemblyfolder/BuildingSite/Contigs/contig_".${$batch_targets_ref}[$i]."/contigs.fa";
                        if (-e $file2) {
                                $cmd .= " ".$file2;
				$c++;
				if (-s $file2 == 0) {
                                        print STDERR "[ERROR] ".$file2." is empty\n";
                                }
                        }

			my $file3 = "$assemblyfolder/BuildingSite/superlocas/superblock_".${$batch_targets_ref}[$i]."/contigs.fa";
			if (-e $file3) {
				$cmd .= " ".$file3;
				$c++;
				if (-s $file3 == 0) {
                                        print STDERR "[ERROR] ".$file3." is empty\n";
                                }
			}
		
			if ($c == 0) {
				print STDERR "[ERROR] Could not find a file for superblock ".${$batch_targets_ref}[$i]." in ".$assemblyfolder.". Was looking for:\n";
				print STDERR "   ", $file1, "\n   ", $file2, "\n   ", $file3, "\n";
			}

		}
	}
	$cmd .= " $leftoverfile > $folder/contigs.fa.tmp";
	system($cmd);

	### Split contigs at stretches of 'N' (fix for problematic Velvet scaffolding)
	### AND
	### Set the fasta header of the left over assemblies to the correct style
	open CONTIGTMP, "$folder/contigs.fa.tmp" or die "Cannot open infile $folder/contigs.fa.tmp\n";
	open CONTIG, ">$folder/contigs.fa" or die "Cannot open outfile $folder/contigs.fa\n";

	my $counter = 1;
	my $lo_counter = 1;
	my $seq = "";
	my $id = "";
	while(<CONTIGTMP>) {
		chomp($_);

		if (substr($_, 0, 1) eq ">") {
			if ($seq ne "") {
				if( ($seq =~ /NNN/) && ($id =~ /velvet/) ) {
					my @sub_seqs = split(/N+/, $seq);
					foreach my $sub_seq (@sub_seqs) {
						print CONTIG ">$counter$id\n$sub_seq\n";
						$counter++;
					}
				}
				else {
					print CONTIG ">$id\n$seq\n";
				}
			}
			if (substr($_, 0, 5) eq ">NODE") {
				$id = check_id(substr($_, 1), $lo_counter);
				$lo_counter++;
			}
			else {
				$id = substr($_, 1);
			}	
			
			$seq = "";
		}
		else {
			$seq .= $_;
		}
	}

	if($seq ne "" && ($seq =~ /NNN/) && ($id =~ /velvet/) ) {
                my @sub_seqs = split(/N+/, $seq);
                foreach my $sub_seq (@sub_seqs) {
                    	print CONTIG ">$counter$id\n$sub_seq\n";
                	$counter++;
	        }
        }
        elsif ($seq ne "") { print CONTIG ">$id\n$seq\n"; }

	close(CONTIGTMP);
	close(CONTIG);

	system("rm $folder/contigs.fa.tmp");

	### Set up reference.fasta file for AMOS
	open REFOUT, ">".$outfolder."/AMOScmp_batches/AMOS_batch_".$batchnum."/ref.fasta";	
	print REFOUT ">".$batchnum."  ".$chr." ".$begin." ".$end."\n";
	print REFOUT substr($REFSEQ{$chr}, $begin-1, ($end-$begin)+1), "\n";
	close REFOUT;

	$batchnum++;
}

sub check_id {
	my ($s, $c) = @_;
	if (substr($s, 0, 4) eq "NODE") {
		return "ID".$c."_velvet_leftover_OZ42";
	}
	return $s;
}


