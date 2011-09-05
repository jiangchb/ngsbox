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
#  Module: Support::Misc::split_error_htm.pl
#  Purpose:
#  In:
#  Out:
#


my $header = "";
my $save_header = "";

# Save the html header #######################################
open(FILE, "Error.htm") || die;
$header = <FILE>;
while( $header !~ />0001</ ) {
	$save_header .= "$header";
	$header = <FILE>;	
}
close(FILE);

# Print html content #########################################
open(FILE, "Error.htm") || die;	
open(OUT, ">Error_header.htm") || die;
while(<FILE>) {
	chomp;

	if( $_ =~ />0001</) {
		print OUT "</table>\n</body>\n</html>";
		open(OUT, ">Error_0001_0050.htm") || die;
		print OUT "$save_header";
	}
	if( $_ =~ />0051</) {
		print OUT "</table>\n</body>\n</html>";
		close( OUT );
		open(OUT, ">Error_0051_0100.htm") || die;
		print OUT "$save_header";
	}
	elsif( $_ =~ />0101</ ) {
		print OUT "</table>\n</body>\n</html>";
		close( OUT );
		open(OUT, ">Error_0101_0150.htm") || die;
		print OUT "$save_header";
	}
	elsif( $_ =~ />0151</ ) {
		print OUT "</table>\n</body>\n</html>";
		close( OUT );
		open(OUT, ">Error_0151_0200.htm") || die;
		print OUT "$save_header";
	}
	
	print OUT "$_\n";
}

open(FILE, "Summary.htm") || die;
open OUTPUT, ">mySummary.htm" or die;
while(<FILE>) {
	if (/<p> <a href='Error.htm' t/) {
		print OUTPUT "<p> <a href='Error_0001_0050.htm'> Error_0001_0050.htm </a> </p>";
		print OUTPUT "<p> <a href='Error_0051_0100.htm'> Error_0051_0100.htm </a> </p>";
		print OUTPUT "<p> <a href='Error_0101_0150.htm'> Error_0101_0150.htm </a> </p>";
		print OUTPUT "<p> <a href='Error_0151_0200.htm'> Error_0151_0200.htm </a> </p>";
		my $skip = <FILE>;
	} else {
		print OUTPUT $_;
	}
}


exit(0);
