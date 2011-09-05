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
#  Module: Analysis::Eighties::genome_matrix::get_full_snp.pl
#  Purpose:
#  In:
#  Out:
#

##########################################################################################
### Read file contain position from command line and get whole 80 ecotype snp/ref call
##########################################################################################

use DBI;


my $usage = "$0 positionfile\n";
my $file = shift or die $usage;
my $dbh;
&connect_to_db();


open FILE, $file or die $usage;

while (my $line = <FILE>) {
        my @a = split " ", $line;


### Read table from database

my $q ="SELECT  *     
        FROM    genome_matrix
        WHERE
                chromosome = $a[0] &&
                position = $a[1] 
#              &&  ref = '$a[2]' 
                  
        ";

my $sth = $dbh->prepare($q);
$sth->execute();


while (my @row = $sth->fetchrow_array) {
                print join("\t", @row);
                print "\n";
        }

#while(my $ref = $sth->fetchrow_hashref()) {
#        print $ref->{position} ."\t"."\n";
#}



}
close(FILE);


#####################################################
### Connects to a database and returns databaseHandle
#####################################################
sub connect_to_db
{
        my $databaseName = "ath_eighty";
        my $driver = "mysql";
        my $host = "orb.eb.local";
        my $username = "jun";
        my $password = "clinnEyt";
        my $dsn = "DBI:$driver:database=$databaseName;host=$host";
        my $drh = DBI->install_driver("mysql");
        $dbh = DBI->connect($dsn, $username, $password ) or die "Could not connect to db. Connect error: $DBI::errstr";
}

