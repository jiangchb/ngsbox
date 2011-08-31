#!/usr/bin/perl
use strict;
use warnings;

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
