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
#  Module: Support::Misc::gc_file_creator.pl
#  Purpose:
#  In:
#  Out:
#

use DBI;

my $dbh;
&connect_to_db();

my $q0 = "SELECT chromosome FROM seq_max ORDER BY chromosome";
my $sth0 = $dbh->prepare($q0);
$sth0->execute();

while(my $ref0 = $sth0->fetchrow_hashref()) {
	my $chr = $ref0->{chromosome};
	print ">$chr\n";

	my $q = "SELECT GC_content_101 FROM seq_ref WHERE chromosome = $chr ORDER BY position";
	my $sth = $dbh->prepare($q);
	$sth->execute();
	while(my $ref = $sth->fetchrow_hashref()) {
		if($ref->{GC_content_101}) {
			print chr(int($ref->{GC_content_101}) + 34);
		}
		else {
			print chr(33);
		}
	}
	$sth->finish();

	print "\n";
}
$sth0->finish();

exit(0);

##############################################################################
### Connects to a database and returns databaseHandle
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
		
