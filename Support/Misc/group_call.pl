#!/usr/bin/perl -Wall
#written by korbinian schneeberger, sept 06

use strict;
use Getopt::Long;

my $working_dir = "";
my $filter = "*";
my %CMD;
my $pretend = 0;

GetCom();

my @dir = glob($working_dir.$filter);

foreach my $f (@dir) {
	my $call = $CMD{call};

	$call =~ s/<tab>/\t/g;
	$call =~ s/<filename>/$f/g;
	$call =~ s/<dirname>/$working_dir/g;
	my @tmp1 = split '/', $f;
	my @tmp = split '\.', $tmp1[$#tmp1];
	$call =~ s/<fn_parsed>/$tmp[0]/g;
	my $fnc = "";
	for (my $i = 0; $i<@tmp-1; $i++) { $fnc .= $tmp[$i]; $fnc .= "." if $i != @tmp-2 }
	$call =~ s/<fn_cut>/$fnc/g;

	if ($pretend == 0) {
  		system($call);
	} else {
		print $call;
	}
}


sub GetCom{

  my @usage = ("Usage: group_call.pl --call=<string> --dir=<dir> --filter=<filter> -p\n

required:
--call\t\tDefine the system call performed for each 
\t\tfile in the directory		                
optional:
--dir\t\tDefine the directory containg the files 
\t\tto be parsed
--filter\tFile filter, default is \"*\"
-p\t\tJust pretend the system call as prints to STDOUT

The call string may contain tokens, which are replaced
by parsed values. Valid tokens are: 
<tab>\t\treplaced by a tab
<filename>\treplaced by the current file name
<dirname>\treplaced by the current dir name
<fn_cut>\treplaced by filename excluding file extension
\t\tBe aware: filename can be empty
<fn_parsed>\treplaced by the 1st part of the 
\t\tfilename til the first dot.\n\n");
 
	die(@usage) if (@ARGV == 0);
	GetOptions(\%CMD, "call=s","dir=s", "p", "filter=s");

	die("Please specify a call string\n") unless $CMD{call};
  
	if (defined($CMD{dir})) {
		$working_dir = $CMD{dir};
		$working_dir .= "/" if substr($working_dir, length($working_dir)-1, 1)  ne "/";
	}
	if (defined($CMD{p})) {
		$pretend = 1;
	}
	if (defined($CMD{filter})) {
		$filter = $CMD{filter} ;
	}
	
	


}
	      

