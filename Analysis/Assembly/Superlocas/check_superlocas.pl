#! /usr/bin/perl

use strict;
my $usage = "$0 superlocasfolder delete\n";
my $folder = shift or die $usage;
my $delete_flag = "";
$delete_flag = shift;
if ($delete_flag eq "") { die $usage; }

chdir($folder);

my @files = glob("superblock_*contigs.fasta");

foreach my $id (@files) { 
	$id =~ s/superblock_//g;
	$id =~ s/contigs.fasta//g;
	my $cmd1 = ("mv superblock_".$id."contigs.fasta superblock_".$id."/contigs.fasta");
	my $cmd2 = ("mv superblock_".$id."contigs._more.fasta superblock_".$id."/contigs._more.fasta");
	my $cmd3 = ("mv superblock_".$id."_pather.info superblock_".$id."/_pather.info");
	my $cmd4 = ("mv superblock_".$id."_path_graph.out superblock_".$id."/_path_graph.out");
	my $cmd5 = ("mv superblock_".$id."_reduced.agraph superblock_".$id."/_reduced.agraph");
	
	system($cmd1);
	system($cmd2);
	system($cmd3);
	system($cmd4);
	system($cmd5);
}


my @folders = glob("superblock_*");

foreach my $folder (@folders) {
	if (-d "$folder") {
		# Catch if contigs.fa was not created correctly:
		my $renew = 0;
		my $fasta_count = 0;
		if (-e $folder."/contigs.fasta") {
			$fasta_count = `grep -c "^>" $folder/contigs.fasta`;
			chomp($fasta_count);
			print STDERR "$folder/contigs.fasta is empty\n" if $fasta_count == 0;
		}
		if (-e $folder."/contigs.fa") {
                        my $c2 = `grep -c "^>" $folder/contigs.fa`;
			chomp($c2);
			if ($fasta_count != $c2) { $renew = 1; }
		}
		# (Re-)make contigs.fa 
		if (!(-e $folder."/contigs.fa") or $renew == 1) {
			if (-e $folder."/contigs.fasta") {
				open FILE, $folder."/contigs.fasta";
				open OUT, ">".$folder."/contigs.fa";
				my $id = "";
				my $count = 0;
				my $seq = "";

				while (my $line = <FILE>) {
					chomp($line);
					if (substr($line, 0, 1) eq ">") {
						if ($seq ne "") {
							$count++;
							print OUT ">ID", $count, "_superlocas_OZ42 | ", $id, "\n", $seq, "\n";
							$seq = "";
						}
						$id = $line;
					}	
					else {
						$seq .= $line;
					}		
				}
				if ($seq ne "") {
                                	$count++;
                                        print OUT ">ID", $count, "_superlocas_OZ42 | ", $id, "\n", $seq, "\n";
                                        $seq = "";
                                }

				close FILE;
				close OUT;
			}
		}
	


		if (!(-e $folder."/log")) {
 	               #print STDERR $folder."/log is missing\n";
                }
		if (!(-e $folder."/_pather.info")) {
                       print STDERR $folder."/_pather.info is missing\n";
                }
		if (!(-e $folder."/contigs.fasta")) {
                	print STDERR $folder."/contigs.fasta is missing\n";
                }

		if ($delete_flag == 0) {	
			# Check existance
			if (!(-e $folder."/contigs._more.fasta")) {
				print STDERR $folder."/contigs._more.fasta is missing\n";
			}
			if (!(-e $folder."/_path_graph.out")) {
				print STDERR $folder."/_path_graph.out is missing\n";
			}
			if (!(-e $folder."/_reduced.agraph")) {
				print STDERR $folder."/_reduced.agraph is missing\n";
			}
		}
		else {
			system("rm -rf ".$folder."/contigs._more.fasta ".$folder."/_path_graph.out ".$folder."/_reduced.agraph");
		}
	}
}




