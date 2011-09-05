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
#  Module: Support::qsub::batch_db_upload.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "\n$0 type<S|I> infolder\n\n";
my $type      = shift or die $usage;
my $infolder  = shift or die $usage;

my @folders = glob($infolder . "/*");

foreach my $folder (@folders) {

	my @folderpath = split("/", $folder);
	my $folderleaf = $folderpath[$#folderpath];

	if( ($type eq "S") && (-e "$folder/SNP_Intersection/merged.all.vcf") ) {

		print "perl /users/GD/tools/pgsp/Support/VCF/vcf2db/vcf2db.pl $folderleaf $folder/SNP_Intersection/merged.all.vcf > $folder/SNP_Intersection/database_upload.txt\n";
		print "sed -i -e 's/SHORE-SHORE/SHORE/g' $folder/SNP_Intersection/database_upload.txt\n";

	}
	elsif( ($type eq "I") && (-e "$folder/Indel_Intersection/merged.union.vcf") ) {

		print "perl /users/GD/tools/pgsp/Support/VCF/vcf2db/vcf2db.pl $folderleaf $folder/Indel_Intersection/merged.union.vcf > $folder/Indel_Intersection/database_upload.txt\n";
		print "sed -i -e 's/SHORE-SHORE/SHORE/g' $folder/Indel_Intersection/database_upload.txt\n";

	}
}
