#!/usr/bin/perl
#written by ks

use strict;
use warnings;
use Getopt::Long;

my %CMD;

GetCom();

my %IDS = ();
while (<FILE2>) {
	my @a = split " ";
	$IDS{$a[0]} = 1;
}


while (<FILE1>) {
        my @a = split " ";
	if (not defined $IDS{$a[0]}) {
		print $_;
	}
}

sub GetCom{

  my @usage = ("\nUsage: drop_lines.pl --file1=<string> --file2=<string>
required:
\t\tfile1	Superset file
\t\tfile2	Set of lines to be dropped in the first file

All lines in file1 not being present in file2 will be parsed out and 
therefore create a distinct set to file2. Lines are defined by their 
first column. Output is written to STDOUT.\n
");
 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "file1=s","file2=s");

	die("Please specify files\n") unless $CMD{file1};
	die("Please specify files\n") unless $CMD{file2};
  
	open FILE1, $CMD{file1} or die "Cannot open ".$CMD{file1}." for reading.\n";	
	open FILE2, $CMD{file2} or die "Cannot open ".$CMD{file2}." for reading.\n";
	
}
	      

