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
#  Module: Support::Database::parse_segments.pl
#  Purpose:
#  In:
#  Out:
#


use Getopt::Long;
use DBI;

my $dbh;

my %CMD = ();
my $database;
my $table;
my $condition;
my $report1;
my $report2;
my $chromosome;
my $position;
my $agglomerator1_method;
my $agglomerator1_attribute;
my $agglomerator2_method;
my $agglomerator2_attribute;
my $agglomerator3_method;
my $agglomerator3_attribute;


GetCom();
&connect_to_db();

my $q2 = "SELECT $chromosome, $position";
if (defined($report1)) {
 	$q2 .= ", $report1";
} 
if (defined($report2)) {
	$q2 .= ", $report2";
}
if (defined($agglomerator1_attribute)) {
	$q2 .= ", $agglomerator1_attribute";
}
if (defined($agglomerator2_attribute)) {
        $q2 .= ", $agglomerator2_attribute";
}
if (defined($agglomerator3_attribute)) {
        $q2 .= ", $agglomerator3_attribute";
}
$q2 .= " FROM $table";
if (defined($condition)) {
	$q2 .= " where $condition";
} 
$q2 .= " order by $chromosome, $position";
print STDERR "Asking db for:\n";
print STDERR "$q2\n";
my $sth2 = $dbh->prepare($q2); $sth2->execute();

my $last_chr = -1;
my $last = -1;
my $start = -1;
my $report_string1 = "";
my $report_string2 = "";

my $agglomerator1_num = 0;
my $agglomerator1_sum = 0;
my $agglomerator1_max = "";
my $agglomerator1_min = "";

my $agglomerator2_num = 0;
my $agglomerator2_sum = 0;
my $agglomerator2_max = "";
my $agglomerator2_min = "";

my $agglomerator3_num = 0;
my $agglomerator3_sum = 0;
my $agglomerator3_max = "";
my $agglomerator3_min = "";

while (my $res = $sth2->fetchrow_hashref()) {
	my $chr = $res->{$chromosome};
	my $pos = $res->{$position};
	my $rep1 = $res->{$report1} if defined($report1);
	my $rep2 = $res->{$report2} if defined($report2);
	my $agg1 = $res->{$agglomerator1_attribute} if defined($agglomerator1_attribute);
	my $agg2 = $res->{$agglomerator2_attribute} if defined($agglomerator2_attribute);
	my $agg3 = $res->{$agglomerator3_attribute} if defined($agglomerator3_attribute);

	if ($pos-1 != $last or $chr != $last_chr) {
		if ($start != -1) {
			print $last_chr, "\t", $start, "\t", $last;
			if (defined($agg1)) {
				print "\t", $agglomerator1_max if ($agglomerator1_method eq "max");
				print "\t", $agglomerator1_min if ($agglomerator1_method eq "min");
				print "\t", $agglomerator1_sum if ($agglomerator1_method eq "sum");
				print "\t", $agglomerator1_sum/$agglomerator1_num if ($agglomerator1_method eq "avg");
			}
			if (defined($agg2)) {
                                print "\t", $agglomerator2_max if ($agglomerator2_method eq "max");
                                print "\t", $agglomerator2_min if ($agglomerator2_method eq "min");
                                print "\t", $agglomerator2_sum if ($agglomerator2_method eq "sum");
                                print "\t", $agglomerator2_sum/$agglomerator2_num if ($agglomerator2_method eq "avg");
                        }
			if (defined($agg3)) {
                                print "\t", $agglomerator3_max if ($agglomerator3_method eq "max");
                                print "\t", $agglomerator3_min if ($agglomerator3_method eq "min");
                                print "\t", $agglomerator3_sum if ($agglomerator3_method eq "sum");
                                print "\t", $agglomerator3_sum/$agglomerator3_num if ($agglomerator3_method eq "avg");
                        }
			if ($report_string1 ne "") {
                                print "\t", $report_string1;
                        }
                        if ($report_string2 ne "") {
                                print "\t", $report_string2;
                        }
			print "\n";
		}

		$start = $pos;
		$report_string1 = "";
		$report_string2 = "";

		$agglomerator1_num = 0;
		$agglomerator1_sum = 0;
		$agglomerator1_max = "";
		$agglomerator1_min = "";

		$agglomerator2_num = 0;
                $agglomerator2_sum = 0;
                $agglomerator2_max = "";
                $agglomerator2_min = "";
	
		$agglomerator3_num = 0;
                $agglomerator3_sum = 0;
                $agglomerator3_max = "";
                $agglomerator3_min = "";
	}


	if (defined($agg1)) {
		$agglomerator1_num++;
		$agglomerator1_sum += $agg1;
		if ($agglomerator1_max eq "" or $agg1 > $agglomerator1_max) {
			$agglomerator1_max = $agg1;
		}
		if ($agglomerator1_min eq "" or $agg1 < $agglomerator1_min) {
			$agglomerator1_min = $agg1;
		}
	}
	if (defined($agg2)) {
                $agglomerator2_num++;
                $agglomerator2_sum += $agg2;
                if ($agglomerator2_max eq "" or $agg2 > $agglomerator2_max) {
                        $agglomerator2_max = $agg2;
                }
                if ($agglomerator2_min eq "" or $agg2 < $agglomerator2_min) {
                        $agglomerator2_min = $agg2;
                }
        }
	if (defined($agg3)) {
                $agglomerator3_num++;
                $agglomerator3_sum += $agg3;
                if ($agglomerator3_max eq "" or $agg3 > $agglomerator3_max) {
                        $agglomerator3_max = $agg3;
                }
                if ($agglomerator3_min eq "" or $agg3 < $agglomerator3_min) {
                        $agglomerator3_min = $agg3;
                }
        }

	$report_string1 .= "," if $report_string1 ne "";
	$report_string1 .= $rep1 if defined($rep1);
	$report_string2 .= "," if $report_string2 ne "";
        $report_string2 .= $rep2 if defined($rep2);
	$last = $pos;
	$last_chr = $chr;
}



exit(0);


sub GetCom {

  my @usage = ("\nUsage: $0

required:
--database\tDB name
--table\t\tDB table name
--chromosome\tAttribute name of db table descibing the chromosomal information
--position\tAttribute name of db table describing the positional information 

optional:
--condition\tUsed as a mysql condition in the WHERE clause (e.g. \"coverage >= 3 and coverage <=25\")

--report1\tReport every value of an attribute in a segment positionwise 
--report2\tDito for a second attribute

--agglomerator1\tSet an agglomeration of a numeric db attribute: Possible values are \"max\", \"min\", \"sum\", \"avg\"
\t\tand add the attribute name in brackets. (e.g. --agglomerator1=\"max(coverage)\")
--agglomerator2\tDito for a second attribute
--agglomerator2\tDito for a third attribute

\n");
        die(@usage) if (@ARGV == 0);
        GetOptions(\%CMD, "database=s", "table=s", "condition=s", "position=s", "chromosome=s", "report1=s", "report2=s", "agglomerator1=s", "agglomerator2=s", "agglomerator3=s");

	die("Please specify database\n") unless defined($CMD{database});
	die("Please specify db table\n") unless defined($CMD{table});
	die("Please specify chromosome\n") unless defined($CMD{chromosome});
	die("Please specify position\n") unless defined($CMD{position});

	$database = $CMD{database};
	$position = $CMD{position};
	$chromosome = $CMD{chromosome};
	$table = $CMD{table};

	if (defined($CMD{condition})) {
	        $condition = $CMD{condition};
	}
	if (defined($CMD{report1})) {
                $report1 = $CMD{report1};
        }
	if (defined($CMD{report2})) {
                $report2 = $CMD{report2};
        }
	if (defined($CMD{agglomerator1})) {
		my @a = split '\(', $CMD{agglomerator1};
		$agglomerator1_method = $a[0];
		$a[1] =~ s/\)//g;
		$agglomerator1_attribute = $a[1];
	}
	if (defined($CMD{agglomerator2})) {
                my @a = split '\(', $CMD{agglomerator2};
                $agglomerator2_method = $a[0];
                $a[1] =~ s/\)//g;
                $agglomerator2_attribute = $a[1];
        }
	if (defined($CMD{agglomerator3})) {
                my @a = split '\(', $CMD{agglomerator3};
                $agglomerator3_method = $a[0];
                $a[1] =~ s/\)//g;
                $agglomerator3_attribute = $a[1];
        }

	return(0);
}

sub connect_to_db {
        my $databaseName = $database;
        my $driver = "mysql";
        my $host = "ume.fml.local";
        my $username = "solexa";
        my $password = "s0lexa";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}


exit(0);
