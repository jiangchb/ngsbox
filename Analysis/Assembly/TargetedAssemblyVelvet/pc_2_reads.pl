#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use lib "$ENV{PGSP}/Assembly";
use parse_mapping;

use lib "$ENV{PGSP}/Prediction/common";
use mask;

my $chrom = shift;
my $flanking_region_length = shift;
my $ecotype = shift;
my $region_file  = shift;
my $mapping_file = shift;
my $outdir = shift;

my $dbh;
&connect_to_db();

open IN, "$region_file" or die "Cannot open pc region file\n";
open FLAT,   ">$outdir/reads_0.fl" or die "Cannot open outfile\n";
open FILTER, ">$outdir/filter.txt" or die "Cannot filter file\n";
open TABLE,  ">$outdir/positions.txt"  or die "Cannot open table file\n";

my $last_chr = -500;
my $last_pos = -500;
my $last_end = -500;

while(<IN>) {
	chomp;
	my ($chr, $pos, $end) = split (/\t/, $_);

	if( ( $pos > ($last_end + (2 * $flanking_region_length)) ) || ( $chr != $last_chr ) ) {
		if ($last_chr == $chrom) {
			&write_pc();
		}
		$last_pos = $pos;
	}

	if( ($end >= $last_end) || ($chr != $last_chr) ) {
		$last_end = $end;
	}
	$last_chr = $chr;
}
if ($last_chr == $chrom) {
	&write_pc();
}



exit(0);


#####################################################################
### Write all reads surrounding a PC to flat file format (reads_0.fl)
sub write_pc
{
	if($last_chr == $chrom) {
		my $mask = new mask($last_chr, $last_pos, $last_end);
		$mask->get_oversampled($dbh, "poly_oversampled_joined", $ecotype, 5);
		$mask->get_repeats($dbh, "repeat_segments", $ecotype, 0);
		my $os_flag = 0; my $repeat_count = 0;
		foreach my $id ( keys %{$mask->{mask}} ) {
			if   ( $mask->{mask}{$id} eq "OS" ) { $os_flag = 1; }
			elsif( $mask->{mask}{$id} eq "RP" ) { $repeat_count++; }
		}


		my $pc_start = $last_pos;
		my $pc_end = $last_end;
		$last_pos -= $flanking_region_length;
		$last_end += $flanking_region_length;

		my $length = $last_end - $last_pos + 1;

		if(! ( ($os_flag == 1) || ( ($repeat_count / $length) > 0.2) ) ) {
			my @results = parse_mapping::get($mapping_file, $last_chr, $last_pos, $last_end - 35);
			foreach ( @results ) {
				my @elem = split(/\t/, $_);
				if( $elem[6] == 1 ) {
					my $seq = "";
					for( my $i = 0; $i < length($elem[2]); $i++ ) {
						if(substr($elem[2], $i, 1) eq "[") {
							$i+=2;
							$seq .= substr($elem[2], $i, 1);
							$i++;
						}
						else {
							$seq .= substr($elem[2], $i, 1);
						}
					}
					$seq =~ s/-//g;
					#print TABLE "$last_chr\t$last_pos\t$last_end\t$elem[0]\t$elem[1]\t$elem[3]-pc\t$pc_start\t$pc_end\n";
					print TABLE $elem[3], "\t", $elem[0], "\t", $elem[1], "\n";
					print FLAT "$elem[3]\t$seq\t$elem[9]\t$elem[10]\t$elem[11]\n";
				}
			}
		}
		else { print FILTER "$last_chr\t$last_pos\t$last_end\n"; }
	}
}


#####################################################
### Connects to a database and returns databaseHandle
#####################################################
sub connect_to_db
{
        my $databaseName = "solexa";
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

