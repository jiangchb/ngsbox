#!/usr/bin/perl
######################################################################################
#Author 	Stephan Ossowski, Korbinian 
#Date 		05/15/07
#Version	0.9
#Input		Fragment file
#Function	Map a set of short sequence reads to a reference geneome allowing gaps
######################################################################################

package mapping_splicesite;

use strict;
use warnings;

######################################################################################
# Create a fasta file from a solexa sequence flatfile
sub flat2fasta
{
	my $stage = shift;
	my $dir = shift;
	open (FRAGS, "$dir/reads_$stage.fl") or die "cannot open $dir/reads_$stage.fl\n";
	open (FASTA, ">$dir/reads_$stage.fa") or die "cannot open $dir/reads_$stage.fa\n";
	while( <FRAGS> ) {
		chomp;
		my ($frag_name, $frag_seq) = split(/\t/, $_);
		print FASTA ">" . $frag_name . "\n" . $frag_seq . "\n";
	}
	close(FASTA);
	close(FRAGS);
	return(0);
}


######################################################################################
# Compare a Vmatch mapping output with original input file
sub comp
{
	my $stage = shift;
	my $dir = shift;
	my $next_stage = $stage + 1;

	### Open map file and read first map line
	open (MAP, "$dir/map_$stage.vm") or die "Zcannot open $dir/map_$stage.vm\n";
	my $map_line = <MAP>;
	if(not defined $map_line) {
		close MAP or die "Cannot close file\n";
		system("cp $dir/reads_$stage.fl $dir/reads_$next_stage.fl");
		return(0);
	}
	chomp($map_line);
	my @map_elem = split(/\t/, $map_line);

	### Open current and next reads file	
	open (FRAGS, "$dir/reads_$stage.fl") or die "cannot open $dir/reads_$stage.fl\n";
        open (NEXT_FRAGS, ">$dir/reads_$next_stage.fl") or die "cannot open $dir/reads_$next_stage.fl\n";
	
	### Loop through original reads and rewrite reads which are not yet mapped
	while( my $frag_line = <FRAGS> ) {
		chomp($frag_line);
		my ($frag_name, $frag_seq) = split(/\t/, $frag_line);
		
		if($map_elem[0] eq $frag_name) {
			while($map_elem[0] eq $frag_name) {
				$map_line = <MAP>;
				if(not defined $map_line) { $map_elem[0] = "endoffile"; last; }
				chomp($map_line);
				@map_elem = split(/\t/, $map_line);
			}
		}
		else {
			print NEXT_FRAGS "$frag_name\t$frag_seq\n";
		}
	}
	
	### Finish
	close FRAGS or die "Cannot close file\n";
	close MAP or die "Cannot close file\n";
	close NEXT_FRAGS or die "Cannot close file\n";
	return(0);
}


######################################################################################
### Parse Vmatch output and write to map format (one row per hit)
sub clean_vmatch_hamming_distance
{
	my $stage          = shift;
	my $dir            = shift;
	my $end_gap        = shift;
	my $repeat_mapping = shift;
	
	open (MAP_ORIGINAL, "$dir/map_original_$stage.vm") or die "cannot open $dir/map_original_$stage.vm\n";
	open (MAP, ">$dir/map_unsorted_$stage.vm") or die "cannot open  $dir/map_unsorted_$stage.vm\n";

	my $junk = <MAP_ORIGINAL>;

	while(my $line =  <MAP_ORIGINAL> ) {
		chomp($line);
		$line =~ s/^\s+//;

		# target_length 0, chr 1, position 2, orientation 3, query_length 4, 
		# frag_name 5, offset 6, vm_mismatches 7, evalue 8
		my @map_elem = split(/\s+/, $line);
		my $target_length = $map_elem[0];
		my $chr = $map_elem[1];
		my $position = $map_elem[2] + 1;
		my $orientation = $map_elem[3];
		my $query_length = $map_elem[4];
		my $frag_name = $map_elem[5];
		my $offset = $map_elem[6];
		my $vm_mismatches = $map_elem[7];
		my $evalue = $map_elem[8];
		$vm_mismatches = abs($vm_mismatches);

		# Read Alignment
		my $alignment = <MAP_ORIGINAL>;
		chomp($alignment);

		if($end_gap > 0) {
			# frag_name, evalue, chr, position, orientation, vm_mismatches,
			# query_length, offset, alignment
			print MAP "$frag_name\t$evalue\t$chr\t$position\t$orientation\t" .
				  "$vm_mismatches\t$query_length\t$offset\t$alignment\n";
		}
		else {
			# frag_name, chr, position, orientation, vm_mismatches,
			# query_length, offset, alignment
			print MAP "$frag_name\t$chr\t$position\t$orientation\t" .
				  "$vm_mismatches\t$query_length\t$offset\t$alignment\n";
		}
	}

	close MAP; close MAP_ORIGINAL;

	if($end_gap > 0) {
		system("sort -g -k1 -k2 --buffer-size=30% $dir/map_unsorted_$stage.vm | uniq > $dir/map_$stage.vm.sorted");
	}
	else {
		system("sort -n --buffer-size=30% $dir/map_unsorted_$stage.vm | uniq > $dir/map_$stage.vm");
	}

	system("rm $dir/map_original_$stage.vm $dir/map_unsorted_$stage.vm");
	
	return(0);
}


######################################################################################
### Parse Vmatch output and write to map format (one row per hit)
sub clean_vmatch_edit_distance
{
	my $stage           = shift;
	my $dir             = shift;
	my $min_read_length = shift;

	open (MAP_ORIGINAL, "$dir/map_original_$stage.vm") or die "cannot open $dir/map_original_$stage.vm\n";
	open (MAP, ">$dir/map_unsorted_$stage.vm") or die "cannot open  $dir/map_unsorted_$stage.vm\n";

	my $junk = <MAP_ORIGINAL>;
	
	### Parse Vmatch output
	while(my $line =  <MAP_ORIGINAL> ) {
		chomp($line);
		$line =~ s/^\s+//;

		# target_length 0, chr 1, position 2, orientation 3, query_length 4, 
		# frag_name 5, offset 6, vm_mismatches 7, evalue 8
		my @map_elem = split(/\s+/, $line);
		my $target_length = $map_elem[0];
		my $chr = $map_elem[1];
		my $position = $map_elem[2] + 1;
		my $orientation = $map_elem[3];
		my $query_length = $map_elem[4];
		my $frag_name = $map_elem[5];
		my $offset = $map_elem[6];
		my $vm_mismatches = $map_elem[7];
		my $evalue = $map_elem[8];
		$vm_mismatches = abs($vm_mismatches);

		if( $query_length < $min_read_length ) {
			<MAP_ORIGINAL>;
			$junk = <MAP_ORIGINAL>;
			if($junk !~ /Query/) { <MAP_ORIGINAL>; }
			next;
		}

		### Read Alignment
		my $alignment = '';
		my $subject = <MAP_ORIGINAL>;
		<MAP_ORIGINAL>;
		my $query = <MAP_ORIGINAL>;
		my ($junk1, $subject_seq, $junk2) = split(/\s/, $subject);
		my ($junk3, $query_seq, $junk4) = split(/\s/, $query);
		for(my $x = 0; $x < length($subject_seq); $x++) {
			my $s = substr($subject_seq, $x, 1);
			my $q = substr($query_seq, $x, 1);
			if($s eq $q) { $alignment .= "$s"; }
			else { $alignment .= "[$s$q]"; }
		}

		# frag_name, evalue, chr, position, orientation, vm_mismatches,
		# query_length, offset, alignment
		print MAP "$frag_name\t$evalue\t$chr\t$position\t$orientation\t" .
			  "$vm_mismatches\t$query_length\t$offset\t$alignment\n";
	}

	close MAP; close MAP_ORIGINAL;

	system("sort -g -k1 -k2 --buffer-size=30% $dir/map_unsorted_$stage.vm | uniq > $dir/map_$stage.vm.sorted");
	system("rm $dir/map_original_$stage.vm $dir/map_unsorted_$stage.vm");

	return(0);
}


######################################################################################
### Drop all reads with higher than min evalue
sub best_evalue
{
	my $stage          = shift;
	my $dir            = shift;

	my $last_frag_name = "NULL";
	my $min_evalue     = 9999;

	open( MAP_IN,  "$dir/map_$stage.vm.sorted") or die "cannot open $dir/map_$stage.vm.sorted\n";
	open( MAP_OUT, ">$dir/map_$stage.vm") or die "Xcannot open $dir/map_$stage.vm\n";

	while( <MAP_IN> ) {
		chomp;
		my @map_elem = split(/\t/, $_);
	
		if($map_elem[0] eq $last_frag_name) {
			if($map_elem[1] eq $min_evalue) {
				print MAP_OUT 	"$map_elem[0]\t$map_elem[2]\t$map_elem[3]\t$map_elem[4]\t$map_elem[5]\t" .
						"$map_elem[6]\t$map_elem[7]\t$map_elem[8]\n";
			}
		}
		else {
			$min_evalue = $map_elem[1];
			print MAP_OUT   "$map_elem[0]\t$map_elem[2]\t$map_elem[3]\t$map_elem[4]\t$map_elem[5]\t" . 
					"$map_elem[6]\t$map_elem[7]\t$map_elem[8]\n";
		}

		$last_frag_name = $map_elem[0];
	}

	close MAP_OUT; close MAP_IN;

	#system("rm $dir/map_$stage.vm.sorted");
	
	return(0);
}


######################################################################################
### Concatenated all map-files and sort final map
sub concatenate
{
	my $dir = shift;
	my $stage = shift;
	my $command = "";
	for(my $i = 0; $i <= $stage; $i++) {
		$command .= "$dir/map_$i.vm ";
	}
	system("cat $command > $dir/map_all.vm");
	system("sort --buffer-size=30% $dir/map_all.vm > $dir/map.vm");
	system("rm $command $dir/map_all.vm $dir/*.fa");
}


######################################################################################
# Count hits against the genome of all reads and write file map.name
sub occurrencies
{
	my $dir = shift;
	my $repeat_mapping  = shift;

	my $last_entry = "";
	my @frag_hits = ();
		
	open (MAP, "$dir/map.vm") or die "Cannot open $dir/map.vm\n";
	open (MAP_LIST, ">$dir/map.list.unsorted") or die "cannot open $dir/map.list.unsorted\n";
	
	while( <MAP> ) {
		chomp;
		my $line = $_;

		# frag_name 0, chr 1, position 2, orientation 3, vm_mismatches 4, 
		# query_length 5, offset 6, align 7, prb 8, qCal 9, chas 10
		my @entries = split(/\t/, $_);

		my $current_entry = $entries[0];
		if($current_entry ne $last_entry) {
			my $hits = scalar(@frag_hits);
			foreach my $frag_hit(@frag_hits) {
				my @this_entries = split(/\t/, $frag_hit);

				if( ($repeat_mapping == 1) || ($hits == 1) ) {
					# chr, pos, alignment, frag_name, orientation, 
					# mismatch, hits, length, offset, prb, qCal, chas, original-seq
					print MAP_LIST 	"$this_entries[1]\t$this_entries[2]\t" .
							"$this_entries[7]\t$this_entries[0]\t" .
							"$this_entries[3]\t$this_entries[4]\t" .
							"$hits\t$this_entries[5]\t" .
							"$this_entries[6]\t$this_entries[8]\t" .
							"$this_entries[9]\t$this_entries[10]\t" .
							"$this_entries[11]\n";
				}
			}
			@frag_hits = ();
		}
		
		push(@frag_hits, $line);
		$last_entry = $current_entry;
	}

	my $hits = scalar(@frag_hits);
	foreach my $frag_hit(@frag_hits) {
		my @this_entries = split(/\t/, $frag_hit);
	
		if( ($repeat_mapping == 1) || ($hits == 1) ) {
			# chr, pos, alignment, frag_name, orientation,
			# mismatch, hits, length, offset, prb
			print MAP_LIST  "$this_entries[1]\t$this_entries[2]\t" .
					"$this_entries[7]\t$this_entries[0]\t" .
					"$this_entries[3]\t$this_entries[4]\t" .
					"$hits\t$this_entries[5]\t" .
					"$this_entries[6]\t$this_entries[8]\t" .
					"$this_entries[9]\t$this_entries[10]\t" .
					"$this_entries[11]\n";
		}	
	}

	close(MAP); close(MAP_LIST);

	system("sort -n -k1 -k2 --buffer-size=30% $dir/map.list.unsorted > $dir/map.list");
	system("rm $dir/map.list.unsorted");
	#system("rm $dir/map.list.unsorted $dir/map.vm");
	
	return(0);
}


######################################################################################
### Add probability values to mapping results
sub add_prb
{
	my $dir = shift;
        open (FRAGS, "$dir/reads_0.fl") or die "cannot open $dir/reads_0.fl\n";
        open (MAP, "$dir/map.vm") or die "cannot open  $dir/map.vm\n";
	open (OUT, ">$dir/tmp") or die "cannot open $dir/tmp\n";

	# Read first map line
        my $map_line = <MAP>;
        chomp($map_line);
        my @map_elem = split(/\t/, $map_line);

	# Loop through original reads and rewrite reads which are not yet mapped
        while( my $frag_line = <FRAGS> ) {
                chomp($frag_line);
                my ($frag_name, $frag_seq, $prb, $qCal, $chas) = split(/\t/, $frag_line);

		if(! defined $prb)  { $prb = "NA"; }
		if(! defined $qCal) { $qCal = "NA"; }
		if(! defined $chas) { $chas = "NA"; }

                if($map_elem[0] eq $frag_name) {
                        while($map_elem[0] eq $frag_name) {
				# frag_name, chr, position, orientation, vm_mismatches, 
				# query_length, offset, evalue, align, prb, qCal, chas, original-seq
				print OUT "$map_line\t$prb\t$qCal\t$chas\t$frag_seq\n";

                                $map_line = <MAP>;
                                if(not defined $map_line) { $map_elem[0] = "endoffile"; last; }
                                chomp($map_line);
                                @map_elem = split(/\t/, $map_line);
                        }
                }
        }

        close FRAGS; close MAP; close OUT;
	system("mv $dir/tmp $dir/map.vm");
        return(0);
}

1;
