#!/usr/bin/perl

use strict;
use warnings;

my $file = shift;
open FILE, $file;

while(<FILE>) {
	my @a = split("\t", $_);

	if($a[0] eq '@SQ') {

		if($a[1] eq "SN:1")  { print "\@SQ\tSN:1\tLN:249250621\tM5:1b22b98cdeb4a9304cb5d48026a85128\n"; }
		if($a[1] eq "SN:2")  { print "\@SQ\tSN:2\tLN:243199373\tM5:a0d9851da00400dec1098a9255ac712e\n"; }
		if($a[1] eq "SN:3")  { print "\@SQ\tSN:3\tLN:198022430\tM5:fdfd811849cc2fadebc929bb925902e5\n"; }
		if($a[1] eq "SN:4")  { print "\@SQ\tSN:4\tLN:191154276\tM5:23dccd106897542ad87d2765d28a19a1\n"; }
		if($a[1] eq "SN:5")  { print "\@SQ\tSN:5\tLN:180915260\tM5:0740173db9ffd264d728f32784845cd7\n"; }
		if($a[1] eq "SN:6")  { print "\@SQ\tSN:6\tLN:171115067\tM5:1d3a93a248d92a729ee764823acbbc6b\n"; }
		if($a[1] eq "SN:7")  { print "\@SQ\tSN:7\tLN:159138663\tM5:618366e953d6aaad97dbe4777c29375e\n"; }
		if($a[1] eq "SN:8")  { print "\@SQ\tSN:8\tLN:146364022\tM5:96f514a9929e410c6651697bded59aec\n"; }
		if($a[1] eq "SN:9")  { print "\@SQ\tSN:9\tLN:141213431\tM5:3e273117f15e0a400f01055d9f393768\n"; }
		if($a[1] eq "SN:10") { print "\@SQ\tSN:10\tLN:135534747\tM5:988c28e000e84c26d552359af1ea2e1d\n"; }
		if($a[1] eq "SN:11") { print "\@SQ\tSN:11\tLN:135006516\tM5:98c59049a2df285c76ffb1c6db8f8b96\n"; }
		if($a[1] eq "SN:12") { print "\@SQ\tSN:12\tLN:133851895\tM5:51851ac0e1a115847ad36449b0015864\n"; }
		if($a[1] eq "SN:13") { print "\@SQ\tSN:13\tLN:115169878\tM5:283f8d7892baa81b510a015719ca7b0b\n"; }
		if($a[1] eq "SN:14") { print "\@SQ\tSN:14\tLN:107349540\tM5:98f3cae32b2a2e9524bc19813927542e\n"; }
		if($a[1] eq "SN:15") { print "\@SQ\tSN:15\tLN:102531392\tM5:e5645a794a8238215b2cd77acb95a078\n"; }
		if($a[1] eq "SN:16") { print "\@SQ\tSN:16\tLN:90354753\tM5:fc9b1a7b42b97a864f56b348b06095e6\n"; }
		if($a[1] eq "SN:17") { print "\@SQ\tSN:17\tLN:81195210\tM5:351f64d4f4f9ddd45b35336ad97aa6de\n"; }
		if($a[1] eq "SN:18") { print "\@SQ\tSN:18\tLN:78077248\tM5:b15d4b2d29dde9d3e4f93d1d0f2cbc9c\n"; }
		if($a[1] eq "SN:19") { print "\@SQ\tSN:19\tLN:59128983\tM5:1aacd71f30db8e561810913e0b72636d\n"; }
		if($a[1] eq "SN:20") { print "\@SQ\tSN:20\tLN:63025520\tM5:0dec9660ec1efaaf33281c0d5ea2560f\n"; }
		if($a[1] eq "SN:21") { print "\@SQ\tSN:21\tLN:48129895\tM5:2979a6085bfe28e3ad6f552f361ed74d\n"; }
		if($a[1] eq "SN:22") { print "\@SQ\tSN:22\tLN:51304566\tM5:a718acaa6135fdca8357d5bfe94211dd\n"; }
		if($a[1] eq "SN:X")  { print "\@SQ\tSN:X\tLN:155270560\tM5:7e0e2e580297b7764e31dbc80c2540dd\n"; }
		if($a[1] eq "SN:Y")  { print "\@SQ\tSN:Y\tLN:59373566\tM5:1fa3474750af0948bdf97d5a0ee52e51\n"; }
		if($a[1] eq "SN:MT") { print "\@SQ\tSN:MT\tLN:16571\tM5:d2ed829b8a1628d16cbeee88e88e39eb\n"; }
	}
	else {
		print $_;
	}
}

close FILE;

exit(0);
