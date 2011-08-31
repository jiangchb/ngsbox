#!/usr/bin/perl
use strict;
use warnings;
use DBI;

my $ecotype    = shift;
my $chr        = shift;
my $min_length = shift;
my $max_length = shift;
my $table      = shift; # poly_pr or poly_unsequenced
my $dbh;
&connect_to_db();

my @chr_length = (0, 30432563, 19705359, 23470805, 18585042, 26992728);


### Load Zeller PRs
my $q ="SELECT 	begin, end 
	FROM 	clark20_prp_ml
	WHERE 	ecotype = '$ecotype' &&
		chromosome = $chr
	ORDER BY chromosome, begin
";
my $sth = $dbh->prepare($q);
$sth->execute();

my %pr_perlegen = ();
while(my $ref = $sth->fetchrow_hashref()) {
	for(my $i = $ref->{begin}; $i <= $ref->{end}; $i++) {
		$pr_perlegen{$i} = 1;
	}
}


### Load Solexa PRs
$q = "	SELECT 	begin, end, length
	FROM	$table
	WHERE 	ecotype = '$ecotype' &&
		chromosome = $chr &&
		length between $min_length AND $max_length
	ORDER by chromosome, begin
";
$sth = $dbh->prepare($q);
$sth->execute();

my %pr_solexa = ();
while(my $ref = $sth->fetchrow_hashref()) {
	my $beg = $ref->{begin};
	my $end = $ref->{end};
	$pr_solexa{$beg} = $end;
}


### Compare
for( my $permute = 1; $permute < 1001; $permute++) {
	my $current_overlap_length = 0;
	my $current_pr_total_length = 0;
	my $curr_rand = int(rand($chr_length[$chr])) + 1;
		
	foreach my $sol_beg (keys %pr_solexa) {
		my $sol_end = $pr_solexa{$sol_beg};
		my $sol_len = $sol_end - $sol_beg + 1;
		$current_pr_total_length += $sol_len;

		my $rand_sol_beg = (($sol_beg + $curr_rand)%$chr_length[$chr]) + 1;
		my $rand_sol_end = (($sol_end + $curr_rand)%$chr_length[$chr]) + 1;

		for(my $pos = $rand_sol_beg; $pos <= $rand_sol_end; $pos++) {
			if(exists $pr_perlegen{$pos}) {
				$current_overlap_length++;
			}
		}
	}

	print "$chr\t$current_overlap_length\t$current_pr_total_length\t" .
		$current_overlap_length/$current_pr_total_length . "\n";
}

exit(0);


#####################################################
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

