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
#  Module: Support::ftp::ftp_get_NCBI.pl
#  Purpose:
#  In:
#  Out:
#

#    ftp_get_NCBI.pl
#    Automatically fetches RefSeq release from NCBI
#    Authors: Joffrey Fitz, Stephan Ossowski
#    Copyright (C) 2005-2009 by Max Planck Institute for Developmental Biology, Tuebingen.
##########################################################################################
use Net::FTP;

local $| = 1;

### Directory structure of refseq/release on ftp.ncbi.nih.gov:
my $host = "ftp.ncbi.nih.gov";
my $dir = "/refseq/release";
my $suffix = ".genomic.fna.gz";
my $tmp_dir = "/ebio/abt6_analysis/nobackup/data/Shore/NCBI_complete/tmp";


### Connect and go to dir
my $ftp = Net::FTP->new($host, Debug => 0) or die "Cannot connect to %host: $@";
$ftp->login("anonymous",'-anonymous@') or die "Cannot login as anonymous", $ftp->message;
$ftp->binary();
$ftp->cwd($dir) or die " ERROR: Cannot change working directory: ", $ftp->message;


### Parse Folders
my @species = $ftp->ls();

foreach my $s (@species) {

	next if($s eq "complete" or $s eq "README");


	### Parse files in folder
	my $subdir = "$dir/$s";
	print "SUBDIR: $subdir\n";
	
	$ftp->cwd($subdir) or die " ERROR: Cannot change working directory: ", $ftp->message;
	my @files = $ftp->ls("*$suffix");

	foreach my $f (@files) {
	
		my $remote_file = "$subdir/$f";
		print "REMOTEFILE: $remote_file\n";
		my $local_file = "$tmp_dir/$subdir/$f";
		print "LOCALFILE: $local_file\n";

		if(stat $local_file) {
			print "  $local_file already exist. Skipping.\n\n";
			next;
		}

		system("mkdir -p $tmp_dir/$subdir");


		### Get file
		if($ftp->get($remote_file, $local_file)) {
			print "$remote_file [OK]\n";
		}
		else {
			print STDERR "[ERROR] Can not stat file: $remote_file\n";
			next;
		}
	}
}

### Finish
$ftp->quit;
exit 0;
