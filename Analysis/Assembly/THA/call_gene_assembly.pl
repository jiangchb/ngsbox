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
#  Module: Analysis::Assembly::THA::call_gene_assembly.pl
#  Purpose:
#  In:
#  Out:
#

~/shore/startshore.sh THA -L libraries.txt -i map.list -l left_over/left_over_missing.all.fl -f /ebio/abt6_analysis/nobackup/data/Shore/Plants/ATH/TAIR8_Masked/TAIR8.v1.Masked.fa.shore -o THA -v -u -y -M 36 -j 23,27,31,35,39,43,47,51 -T /ebio/abt6/korbinian/pgsp/Analysis/Eighties/THA_pipeline/ACD6_Eunyoung.txt -s -e left_over/left_over_missing_paired.lib1.fq,left_over/left_over_single.all.fq

combine contigs von allen assembly tools

