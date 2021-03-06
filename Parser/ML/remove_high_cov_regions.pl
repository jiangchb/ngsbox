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
#  Module: Parser::ML::remove_high_cov_regions.pl
#  Purpose:
#  In:
#  Out:
#


my $usage = "perl remove_high_cov_regions.pl map.list";
my $file = shift or die $usage;
open FILE, $file;

# Nonrep coverage > 100 in Col-0
# and combined regions which are less than 10kb apart

@REG=("1	12294231	12294704",
"1	14476897	14513799",
"1	14547344	14548147",
"1	14584241	14584273",
"1	15042355	15042360",
"1	15084390	15104466",
"1	15207243	15207279",
"1	15439883	15440340",
"1	15589789	15589823",
"1	15692836	15692862",
"1	15749766	15749804",
"1	16434456	16435924",
"1	16515818	16529479",
"2	1035	10372",
"2	2892244	2892251",
"2	2938627	2938653",
"2	3088690	3088699",
"2	3244757	3509725",
"2	3605858	3636522",
"2	3668675	3671389",
"2	3720470	3724808",
"2	4136446	4137591",
"2	6422348	6422381",
"2	7198017	7198373",
"2	17953279	17953441",
"2	19704833	19704868",
"3	2094499	2094536",
"3	12253288	12253318",
"3	12674567	12675468",
"3	13536608	13541390",
"3	13590616	13602399",
"3	13644762	13674824",
"3	13709652	13737666",
"3	13769924	13880922",
"3	14196716	14236124",
"3	15179046	15179254",
"3	15273121	15273122",
"3	15315705	15316365",
"3	15728019	15728558",
"3	15795628	15795642",
"4	2060992	2061028",
"4	2931901	2931938",
"4	2957269	2977003",
"4	3025311	3061412",
"4	3221113	3221142",
"4	3842840	3842856",
"4	3870043	3872504",
"4	3950565	4011497",
"4	4761723	4778405",
"4	11099787	11099832",
"5	3253093	3253238",
"5	4508390	4508425",
"5	11181195	11188898",
"5	11341811	11341847",
"5	11668511	11846350",
"5	11890589	11993038",
"5	12042236	12082446",
"5	12167539	12167573",
"5	12903623	12903990",
"5	15687334	15690785");

my %IGNORE = ();
foreach my $reg(@REG) {
	my @a = split " ", $reg;
	for (my $i=$a[1]; $i<$a[2]-36; $i++) {
		$IGNORE{$a[0]."#".$i} = 1;
	}
}

while (my $line = <FILE>) {
	my @a = split " ", $line;
	if (!defined($IGNORE{$a[0]."#".$a[1]})) {
		print $line;
	}
}


