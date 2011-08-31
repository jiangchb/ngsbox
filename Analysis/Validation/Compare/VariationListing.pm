#!/usr/bin/perl 

use strict;
use warnings;

package VariationListing;

sub new {
	my $self = {
		snp => {},
		del => {},
		ins => {},
		snp_num => 0,
		del_num => 0,
		ins_num => 0,
		snp_num_priv => 0,
                del_num_priv => 0,
                ins_num_priv => 0,
		del_num_priv_1bp => 0,
                del_num_priv_4bp => 0,
                ins_num_priv_1bp => 0,
                ins_num_priv_4bp => 0,
		snp_num_id => 0,
		del_num_id => 0,
                ins_num_id => 0,
		del_num_id_1bp => 0,
		del_num_id_4bp => 0,
		ins_num_id_1bp => 0,
                ins_num_id_4bp => 0,
		snp_overlap => {},
	};
	bless $self;
	return $self;
}

#############################################################################################
### Calculates the genome sequence of a sampled genome
sub calc_private_variation_sets {
	my ($self, $snp1file, $del1file, $ins1file, $snp2file, $del2file, $ins2file, $debug) = @_;

	$self->{del_num_1bp} = 0;
	$self->{del_num_4bp} = 0;	
	$self->{del_num_priv_1bp} = 0;
	$self->{del_num_priv_4bp} = 0;
	$self->{del_num_id_1bp} = 0;
	$self->{del_num_id_4bp} = 0;

        $self->{ins_num_1bp} = 0;
        $self->{ins_num_4bp} = 0;
        $self->{ins_num_priv_1bp} = 0;
        $self->{ins_num_priv_4bp} = 0;
        $self->{ins_num_id_1bp} = 0;
        $self->{ins_num_id_4bp} = 0;


	# SNP file
	open SNP1, $snp1file or die "Cannot open file\n";
	while (my $line = <SNP1>) {
        	my @a = split " ", $line;
	        $self->{snp}{$a[1]."#".$a[2]} = $a[4];
        	$self->{snp_num}++;
	}
	close SNP1;
	$self->{snp_num_priv} = $self->{snp_num};
	open SNP2, $snp2file or die "Cannot open file\n";
	while (my $line = <SNP2>) {
        	my @a = split " ", $line;
	        if (defined($self->{snp}{$a[1]."#".$a[2]})) {
			$self->{snp_overlap}{$a[1]."#".$a[2]} = $a[4];
        	        delete $self->{snp}{$a[1]."#".$a[2]};
                	$self->{snp_num_priv}--;
	        }
	}
	$self->{snp_num_id} = $self->{snp_num} - $self->{snp_num_priv};
	close SNP2;
	#print STDERR "To check: SNPs: $num_snps\n";

	# Deletion
	open DEL1, $del1file or die "Cannot open file\n";
        while (my $line = <DEL1>) {
                my @a = split " ", $line;
                $self->{del}{$a[1]."#".$a[2]."#".$a[3]} = $a[4];
                $self->{del_num}++;
		if ($a[4] <= 3) {
			$self->{del_num_1bp}++;
		}
		else {
			$self->{del_num_4bp}++;
		}
        }
        close DEL1;
	$self->{del_num_priv} = $self->{del_num};
	$self->{del_num_priv_1bp} = $self->{del_num_1bp};
	$self->{del_num_priv_4bp} = $self->{del_num_4bp};	
#print "del: ", $self->{del_num_priv}, " ", $self->{del_num_priv_1bp}, " ", $self->{del_num_priv_4bp} , "\n";
        open DEL2, $del2file or die "Cannot open file\n";
        while (my $line = <DEL2>) {
                my @a = split " ", $line;
                if (defined($self->{del}{$a[1]."#".$a[2]."#".$a[3]})) {
                        delete $self->{del}{$a[1]."#".$a[2]."#".$a[3]};
                        $self->{del_num_priv}--;
			if ($a[4] <= 3) {
				$self->{del_num_priv_1bp}--;
			}
			else {
				$self->{del_num_priv_4bp}--;
			}
                }
        }
	$self->{del_num_id} = $self->{del_num} - $self->{del_num_priv};
	$self->{del_num_id_1bp} = $self->{del_num_1bp} - $self->{del_num_priv_1bp};
	$self->{del_num_id_4bp} = $self->{del_num_4bp} - $self->{del_num_priv_4bp};
#print "del: ", $self->{del_num_id}, " ", $self->{del_num_id_1bp}, " ", $self->{del_num_id_4bp} , "\n";
        close DEL2;
        #print STDERR "To check: DELs: $num_del\n";

	# Insertion
        open INS1, $ins1file or die "Cannot open file\n";
        while (my $line = <INS1>) {
                my @a = split " ", $line;
                $self->{ins}{$a[1]."#".$a[2]."#".$a[3]} = $a[5];
                $self->{ins_num}++;
		if (length($a[5]) <= 3) {
			$self->{ins_num_1bp}++;
		}
		else {
			$self->{ins_num_4bp}++;
		}
        }
        close INS1;
	$self->{ins_num_priv} = $self->{ins_num};
	$self->{ins_num_priv_1bp} = $self->{ins_num_1bp};
	$self->{ins_num_priv_4bp} = $self->{ins_num_4bp};
        open INS2, $ins2file or die "Cannot open file\n";
        while (my $line = <INS2>) {
                my @a = split " ", $line;
                if (defined($self->{ins}{$a[1]."#".$a[2]."#".$a[3]})) {
			if (length($self->{ins}{$a[1]."#".$a[2]."#".$a[3]}) <= 3) {
				$self->{ins_num_priv_1bp}--;
			}
			else {
				$self->{ins_num_priv_4bp}--;
			}
			delete $self->{ins}{$a[1]."#".$a[2]."#".$a[3]};
                        $self->{ins_num_priv}--;
                }
        }
	$self->{ins_num_id} = $self->{ins_num} - $self->{ins_num_priv};
	$self->{ins_num_id_1bp} = $self->{ins_num_1bp} - $self->{ins_num_priv_1bp};
	$self->{ins_num_id_4bp} = $self->{ins_num_4bp} - $self->{ins_num_priv_4bp};
        close INS2;
        #print STDERR "To check: INSs: $num_ins\n";

}


1;



