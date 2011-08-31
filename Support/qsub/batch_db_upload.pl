#! /usr/bin/perl
use strict;
use warnings;

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
