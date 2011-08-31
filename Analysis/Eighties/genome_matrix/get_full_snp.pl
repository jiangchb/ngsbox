#!/usr/bin/perl
##########################################################################################
### Read file contain position from command line and get whole 80 ecotype snp/ref call
##########################################################################################

use strict;
use warnings;
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

