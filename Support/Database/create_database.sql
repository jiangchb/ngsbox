create database your_assay_here;

--- Tables describing reference sequence ---
CREATE TABLE seq_ref (
  chromosome int(1) NOT NULL default '0',
  position int(8) NOT NULL default '0',
  base char(1) default NULL,
  GC_content double default NULL,
  PRIMARY KEY  (chromosome,position)
);

CREATE TABLE seq_ambiguous (
  chromosome int(7) default NULL,
  position int(9) default NULL,
  upac char(1) default NULL,
  PRIMARY KEY  (chromosome, position)
);

CREATE TABLE seq_max (
  assembly varchar(50) NOT NULL default '',
  chromosome int(11) NOT NULL default '0',
  genome char(1) default NULL,
  max_pos int(11) default NULL,
  PRIMARY KEY  (assembly,chromosome)
);


--- Tables for repeat analysis ---
CREATE TABLE repeat_org_36 (
  chromosome int(11) NOT NULL default '0',
  position int(11) NOT NULL default '0',
  PRIMARY KEY  (chromosome,position)
);

CREATE TABLE repeat_scale_0_36 (
  chr int(11) NOT NULL default '0',
  pos int(11) NOT NULL default '0',
  base char(1) default NULL,
  coverage int(11) default NULL,
  count_A int(11) default NULL,
  count_C int(11) default NULL,
  count_G int(11) default NULL,
  count_T int(11) default NULL,
  count_D int(11) default NULL,
  count_N int(11) default NULL,
  average_hits double default NULL,
  averages_mismatches double default NULL,
  mm_0 int(11) default NULL,
  mm_1 int(11) default NULL,
  mm_2 int(11) default NULL,
  mm_3 int(11) default NULL,
  mm_x int(11) default NULL,
  PRIMARY KEY  (chr,pos),
  KEY covidx (coverage)
);

CREATE TABLE repeat_seg_36 (
  chromosome smallint(5) NOT NULL default '0',
  position int(9) NOT NULL default '0',
  base char(1) default NULL,
  coverage int(9) unsigned default NULL,
  count_A int(9) unsigned default NULL,
  count_C int(9) unsigned default NULL,
  count_G int(9) unsigned default NULL,
  count_T int(9) unsigned default NULL,
  count_D int(9) unsigned default NULL,
  count_N int(9) unsigned default NULL,
  average_hits double default NULL,
  average_mismatches decimal(4,3) default NULL,
  mm_0 int(9) unsigned default NULL,
  mm_1 int(9) unsigned default NULL,
  mm_2 int(9) unsigned default NULL,
  mm_3 int(9) unsigned default NULL,
  mm_x int(9) unsigned default NULL,
  core_base char(1) default NULL,
  core_coverage int(9) unsigned default NULL,
  core_count_A int(9) unsigned default NULL,
  core_count_C int(9) unsigned default NULL,
  core_count_G int(9) unsigned default NULL,
  core_count_T int(9) unsigned default NULL,
  core_count_D int(9) unsigned default NULL,
  core_count_N int(9) unsigned default NULL,
  core_average_hits double default NULL,
  core_average_mismatches decimal(4,3) default NULL,
  core_mm_0 int(9) unsigned default NULL,
  core_mm_1 int(9) unsigned default NULL,
  core_mm_2 int(9) unsigned default NULL,
  core_mm_3 int(9) unsigned default NULL,
  core_mm_x int(9) unsigned default NULL,
  nonrep_base char(1) default NULL,
  nonrep_coverage int(9) unsigned default NULL,
  nonrep_count_A int(9) unsigned default NULL,
  nonrep_count_C int(9) unsigned default NULL,
  nonrep_count_G int(9) unsigned default NULL,
  nonrep_count_T int(9) unsigned default NULL,
  nonrep_count_D int(9) unsigned default NULL,
  nonrep_count_N int(9) unsigned default NULL,
  nonrep_average_mismatches decimal(4,3) default NULL,
  nonrep_mm_0 int(9) unsigned default NULL,
  nonrep_mm_1 int(9) unsigned default NULL,
  nonrep_mm_2 int(9) unsigned default NULL,
  nonrep_mm_3 int(9) unsigned default NULL,
  nonrep_mm_x int(9) unsigned default NULL,
  core_nonrep_base char(1) default NULL,
  core_nonrep_coverage int(9) unsigned default NULL,
  core_nonrep_count_A int(9) unsigned default NULL,
  core_nonrep_count_C int(9) unsigned default NULL,
  core_nonrep_count_G int(9) unsigned default NULL,
  core_nonrep_count_T int(9) unsigned default NULL,
  core_nonrep_count_D int(9) unsigned default NULL,
  core_nonrep_count_N int(9) unsigned default NULL,
  core_nonrep_average_mismatches decimal(4,3) default NULL,
  core_nonrep_mm_0 int(9) unsigned default NULL,
  core_nonrep_mm_1 int(9) unsigned default NULL,
  core_nonrep_mm_2 int(9) unsigned default NULL,
  core_nonrep_mm_3 int(9) unsigned default NULL,
  core_nonrep_mm_x int(9) unsigned default NULL,
  PRIMARY KEY  (chromosome,position)
)

--- Table for sRNA consensus ---
CREATE TABLE consensus_001_0029_2 (
  chromosome smallint(5) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  length tinyint(2) unsigned NOT NULL default '0',
  base char(1) default NULL,
  coverage int(9) unsigned default NULL,
  coverage_forward int(9) unsigned default NULL,
  coverage_reverse int(9) unsigned default NULL,
  average_hits double default NULL,
  average_mismatches decimal(4,3) default NULL,
  nonrep_base char(1) default NULL,
  nonrep_coverage int(9) unsigned default NULL,
  nonrep_coverage_forward int(9) unsigned default NULL,
  nonrep_coverage_reverse int(9) unsigned default NULL,
  nonrep_average_mismatches decimal(4,3) default NULL,
  PRIMARY KEY  (chromosome,position,length),
  KEY chr_pos_len_cov_index (chromosome, position, length, coverage),
  KEY length_index (length)
);

--- Table for sRNA consensus new ---
CREATE TABLE consensus_001_0029_2 (
  chromosome smallint(5) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  length tinyint(2) unsigned NOT NULL default '0',
  average_hits double default NULL,
  coverage int(9) unsigned default NULL,
  coverage_fwd int(9) unsigned default NULL,
  coverage_rev int(9) unsigned default NULL,
  nonrep_coverage int(9) unsigned default NULL,
  nonrep_coverage_fwd int(9) unsigned default NULL,
  nonrep_coverage_rev int(9) unsigned default NULL,
  PRIMARY KEY  (chromosome,position,length),
  KEY chr_pos_len_fwd_index (chromosome, position, length, coverage_fwd),
  KEY chr_pos_len_rev_index (chromosome, position, length, coverage_rev),
  KEY length_index (length)
);

--- Table for sRNA segments ---
CREATE TABLE segment_consensus_001_0029_2_fwd (
  length tinyint(2) unsigned NOT NULL default '0',
  orientation char(1) NOT NULL, 
  chromosome smallint(5) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  segment_length smallint(5) unsigned default NULL,
  read_count int(7) unsigned default NULL,
  avg_coverage double default NULL,
  avg_avg_hits double default NULL,
  PRIMARY KEY (length, orientation, chromosome, begin),
  KEY len_ori_chr_beg_end_index (length, orientation, chromosome, begin, end),
  KEY length_index (length)
);


--- Table for transcriptome consensus ---
CREATE TABLE consensus_001_0040_1 (
  chromosome smallint(5) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  length tinyint(2) unsigned NOT NULL default '0',
  base char(1) default NULL,
  coverage int(9) unsigned default NULL,
  coverage_forward int(9) unsigned default NULL,
  coverage_reverse int(9) unsigned default NULL,
  average_hits double default NULL,
  average_mismatches decimal(4,3) default NULL,
  nonrep_base char(1) default NULL,
  nonrep_coverage int(9) unsigned default NULL,
  nonrep_coverage_forward int(9) unsigned default NULL,
  nonrep_coverage_reverse int(9) unsigned default NULL,
  nonrep_average_mismatches decimal(4,3) default NULL,
  PRIMARY KEY  (chromosome,position),
  KEY chr_pos_len_cov_index (chromosome, position, coverage)
);

--- Table for transcriptome segmentation
CREATE TABLE segmentation (
  ecotype varchar(20) NOT NULL default '',
  chromosome smallint(5) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  length int(9) unsigned default NULL,
  avg_cov double default NULL,
  max_cov double default NULL,
  PRIMARY KEY  (ecotype, chromosome, begin, end)
);


--- Mapping summary tables ---
CREATE TABLE consensus_RS106 (
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',

  base char(1) default NULL,
  coverage int(6) unsigned default NULL,
  count_A int(6) unsigned default NULL,
  count_C int(6) unsigned default NULL,
  count_G int(6) unsigned default NULL,
  count_T int(6) unsigned default NULL,
  count_D int(6) unsigned default NULL,
  count_N int(6) unsigned default NULL,
  average_hits double default NULL,
  average_mismatches decimal(4,3) default NULL,
  mm_0 int(6) unsigned default NULL,
  mm_1 int(6) unsigned default NULL,
  mm_2 int(6) unsigned default NULL,
  mm_3 int(6) unsigned default NULL,
  mm_x int(6) unsigned default NULL,

  core_base char(1) default NULL,
  core_coverage int(6) unsigned default NULL,
  core_count_A int(6) unsigned default NULL,
  core_count_C int(6) unsigned default NULL,
  core_count_G int(6) unsigned default NULL,
  core_count_T int(6) unsigned default NULL,
  core_count_D int(6) unsigned default NULL,
  core_count_N int(6) unsigned default NULL,
  core_average_hits double default NULL,
  core_average_mismatches decimal(4,3) default NULL,
  core_mm_0 int(6) unsigned default NULL,
  core_mm_1 int(6) unsigned default NULL,
  core_mm_2 int(6) unsigned default NULL,
  core_mm_3 int(6) unsigned default NULL,
  core_mm_x int(6) unsigned default NULL,

  nonrep_base char(1) default NULL,
  nonrep_coverage int(6) unsigned default NULL,
  nonrep_count_A int(6) unsigned default NULL,
  nonrep_count_C int(6) unsigned default NULL,
  nonrep_count_G int(6) unsigned default NULL,
  nonrep_count_T int(6) unsigned default NULL,
  nonrep_count_D int(6) unsigned default NULL,
  nonrep_count_N int(6) unsigned default NULL,
  nonrep_average_mismatches decimal(4,3) default NULL,
  nonrep_mm_0 int(6) unsigned default NULL,
  nonrep_mm_1 int(6) unsigned default NULL,
  nonrep_mm_2 int(6) unsigned default NULL,
  nonrep_mm_3 int(6) unsigned default NULL,
  nonrep_mm_x int(6) unsigned default NULL,

  core_nonrep_base char(1) default NULL,
  core_nonrep_coverage int(6) unsigned default NULL,
  core_nonrep_count_A int(6) unsigned default NULL,
  core_nonrep_count_C int(6) unsigned default NULL,
  core_nonrep_count_G int(6) unsigned default NULL,
  core_nonrep_count_T int(6) unsigned default NULL,
  core_nonrep_count_D int(6) unsigned default NULL,
  core_nonrep_count_N int(6) unsigned default NULL,
  core_nonrep_average_mismatches decimal(4,3) default NULL,
  core_nonrep_mm_0 int(6) unsigned default NULL,
  core_nonrep_mm_1 int(6) unsigned default NULL,
  core_nonrep_mm_2 int(6) unsigned default NULL,
  core_nonrep_mm_3 int(6) unsigned default NULL,
  core_nonrep_mm_x int(6) unsigned default NULL,

  ref_base char(1) default NULL,
  gc_cont tinyint(3) unsigned default NULL,
  gc_flag tinyint(1) unsigned default NULL,
  organelle_flag tinyint(4) unsigned default NULL,
  exp_cov double default NULL,
  exp_cov_0 double default NULL,
  seg_flag char(1),
  oversampled char(1),

  PRIMARY KEY  (chromosome, position),
  KEY average_hits_index (average_hits),
  KEY coverage_index (coverage),
  KEY nonrep_coverage_index (nonrep_coverage)
);



CREATE TABLE insertion_bur (
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',

  insertion_length tinyint(1) unsigned NOT NULL default '0',
  insertion_position tinyint(1) unsigned NOT NULL default '0',

  support  int(6) unsigned default NULL,
  coverage int(6) unsigned default NULL,
  count_A  int(6) unsigned default NULL,
  count_C  int(6) unsigned default NULL,
  count_G  int(6) unsigned default NULL,
  count_T  int(6) unsigned default NULL,
  count_N  int(6) unsigned default NULL,

  core_support  int(6) unsigned default NULL,
  core_coverage int(6) unsigned default NULL,
  core_count_A  int(6) unsigned default NULL,
  core_count_C  int(6) unsigned default NULL,
  core_count_G  int(6) unsigned default NULL,
  core_count_T  int(6) unsigned default NULL,
  core_count_N  int(6) unsigned default NULL,

  nonrep_support  int(6) unsigned default NULL,
  nonrep_coverage int(6) unsigned default NULL,
  nonrep_count_A  int(6) unsigned default NULL,
  nonrep_count_C  int(6) unsigned default NULL,
  nonrep_count_G  int(6) unsigned default NULL,
  nonrep_count_T  int(6) unsigned default NULL,
  nonrep_count_N  int(6) unsigned default NULL,

  core_nonrep_support  int(6) unsigned default NULL,
  core_nonrep_coverage int(6) unsigned default NULL,
  core_nonrep_count_A  int(6) unsigned default NULL,
  core_nonrep_count_C  int(6) unsigned default NULL,
  core_nonrep_count_G  int(6) unsigned default NULL,
  core_nonrep_count_T  int(6) unsigned default NULL,
  core_nonrep_count_N  int(6) unsigned default NULL,

  PRIMARY KEY  (chromosome, position, insertion_length, insertion_position),
  KEY support_key(support),
  KEY coverage_index (coverage),
  KEY core_support_index(core_support),
  KEY core_coverage_index (core_coverage),
  KEY nonrep_support_index(nonrep_support),
  KEY nonrep_coverage_index (nonrep_coverage),
  KEY core_nonrep_support_index(core_nonrep_support),
  KEY core_nonrep_coverage_index (core_nonrep_coverage)
);


--- Polymorphism tables ---
CREATE TABLE poly_snp (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  reference char(1) default NULL,
  basecall char(1) default NULL,
  based_on varchar(20) default NULL,
  support int(6) unsigned default NULL,
  concordance double unsigned default NULL,
  qval tinyint(2) unsigned default NULL,
  sval tinyint(2) unsigned default NULL,
  cval tinyint(2) unsigned default NULL,
  repeat_type varchar(20) default NULL,
  stype varchar(20) default NULL,

  PRIMARY KEY (chromosome, position, ecotype),
  KEY chr_pos_index(chromosome, position),
  KEY ref_mb_index(reference, basecall),
  KEY ecotype_index (ecotype)
);


CREATE TABLE poly_snp_coding (
  ecotype varchar(20) NOT NULL default '',
  tair_id varchar(11) NOT NULL default '',
  isoform tinyint(1) NOT NULL default 0,
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  gene_pos int(6) unsigned default NULL,
  coding_pos int(6) unsigned default NULL,
  codon_pos tinyint(1) unsigned default NULL,
  ref_base char(1) default NULL,
  sub_base char(1) default NULL,
  syn_nonsyn varchar(10) default NULL,
  ref_aa char(1) default NULL,
  sub_aa char(1) default NULL,

  PRIMARY KEY (ecotype, tair_id, isoform, chromosome, position),
  KEY eco_chr_pos_index(ecotype, chromosome, position),
  KEY chr_pos_index(chromosome, position),
  KEY ecotype_index (ecotype),
  KEY tair_id_index (tair_id),
  KEY isoform_index (isoform),
  KEY syn_nonsyn_index (syn_nonsyn),
  KEY codon_pos_index (codon_pos)
);
  
CREATE TABLE poly_snp_segment (
  ecotype varchar(20) NOT NULL default '',
  tair_id varchar(11) NOT NULL default '',
  isoform tinyint(1) NOT NULL default 0,
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  ref_base char(1) default NULL,
  sub_base char(1) default NULL,
  stype varchar(20) default NULL,

  PRIMARY KEY (ecotype, tair_id, isoform, chromosome, position),
  KEY eco_chr_pos_index(ecotype, chromosome, position),
  KEY chr_pos_index(chromosome, position),
  KEY ecotype_index (ecotype),
  KEY tair_id_index (tair_id),
  KEY isoform_index (isoform),
  KEY seg_type_index (segment_type)
);

CREATE TABLE poly_indel (
  ecotype varchar(20) NOT NULL default '',
  variation  varchar(10) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  begin INT(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  length tinyint(1) unsigned default NULL,
  seq char(3) default NULL,
  based_on varchar(20) default NULL,
  support int(6) unsigned default NULL,
  concordance double unsigned default NULL,
  repeat_type varchar(20) default NULL,
  stype varchar(20) default NULL,

  PRIMARY KEY  (ecotype, variation, chromosome, begin, end),
  KEY chr_pos_index(chromosome, begin),
  KEY ecotype_index (ecotype)
);


CREATE TABLE poly_indel_segment (
  ecotype varchar(20) NOT NULL default '',
  variation  varchar(10) NOT NULL default '',
  tair_id varchar(11) NOT NULL default '',
  isoform tinyint(1) NOT NULL default 0,
  chromosome int(7) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  length tinyint(1) unsigned default NULL,
  seq char(3) default NULL,
  concordance double unsigned default NULL,
  based_on varchar(20) default NULL,
  support int(6) unsigned default NULL,
  repeat_type varchar(20) default NULL,
  stype varchar(12) default NULL,

  PRIMARY KEY (ecotype, tair_id, isoform, variation, chromosome, begin, end),
  KEY eco_chr_pos_index(ecotype, variation, chromosome, begin, end),
  KEY chr_pos_index(chromosome, begin, end),
  KEY ecotype_index (ecotype),
  KEY tair_id_index (tair_id),
  KEY isoform_index (isoform),
  KEY length_index (length),
  KEY seg_type_index (segment_type)
);


CREATE TABLE poly_pr (
  ecotype varchar(30) default NULL,
  chromosome int(7) unsigned default NULL,
  begin int(9) unsigned default NULL,
  end int(9) unsigned default NULL,
  length int(9) unsigned default NULL,
  ambiguous_count int(9) unsigned default NULL,
  gc_avg double default NULL,
  gc_max smallint(6) default NULL,
  gc_flag_sum int(9) default NULL,
  u_count int(9) default NULL,
  m_count int(9) default NULL,
  r_count int(9) default NULL,
  mm_left char(1) default NULL,
  mm_right char(1) default NULL,
  PRIMARY KEY  (ecotype, chromosome, begin, end),
  KEY chr_pos_index(chromosome, begin, end),
  KEY eco_length_index (ecotype, length),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_unsequenced (
  ecotype varchar(30) default NULL,
  chromosome int(7) unsigned default NULL,
  begin int(9) unsigned default NULL,
  end int(9) unsigned default NULL,
  length int(9) unsigned default NULL,
  ambiguous_count int(9) unsigned default NULL,
  gc_avg double default NULL,
  gc_max smallint(6) default NULL,
  gc_flag_sum int(9) default NULL,
  join_flag tinyint(4) default NULL,
  u_count int(9) default NULL,
  m_count int(9) default NULL,
  r_count int(9) default NULL,
  mm_left char(1) default NULL,
  mm_right char(1) default NULL,
  PRIMARY KEY  (ecotype, chromosome, begin, end),
  KEY chr_pos_index(chromosome, begin, end),
  KEY eco_length_index (ecotype, length),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_oversampled (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  end int(9) unsigned default NULL,
  length int(9) unsigned default NULL,
  factor double default NULL,
  PRIMARY KEY  (ecotype,chromosome,position)
);

CREATE TABLE poly_pc (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) NOT NULL default '0',
  end int(9) NOT NULL default '0',
  PRIMARY KEY  (chromosome,position,end,ecotype),
  KEY chr_pos_end_index (chromosome,position, end)
);

CREATE TABLE poly_contig_range (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  PRIMARY KEY (chromosome,begin,end,ecotype),
  KEY chr_begin_end_index (chromosome,begin,end),
  KEY ecotype_index (ecotype)
);


CREATE TABLE poly_consensus_XXX (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  reference char(1) default NULL,
  base_call char(1) default NULL,
  based_on varchar(20) default NULL,
  support int(6) unsigned default NULL,
  used_coverage int(6) unsigned default NULL,
  repeat_type varchar(20) default NULL,
  PRIMARY KEY  (chromosome, position, ecotype),
  KEY chr_pos_index (chromosome,position),
  KEY ref_mb_index (reference),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_reference_Bur0 (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  basecall char(1) default NULL,
  based_on varchar(20) default NULL,
  support int(6) unsigned default NULL,
  concordance double unsigned default NULL,
  qval tinyint(2) unsigned default NULL,
  sval tinyint(2) unsigned default NULL,
  cval tinyint(2) unsigned default NULL,
  repeat_type varchar(20) default NULL,
  stype varchar(20) default NULL,

  PRIMARY KEY (chromosome, position, ecotype),
  KEY chr_pos_index(chromosome, position),
  KEY ref_mb_index(reference, basecall),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_contig(
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  length int(9) unsigned default NULL,
  seq text default NULL,
  PRIMARY KEY (ecotype, chromosome, begin, end),
  KEY chr_pos_end_index(chromosome, begin, end),
  KEY length_index (length)
);

CREATE TABLE poly_allele_snp (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  reference char(1) default NULL,
  allele_1_base char(1) default NULL,
  allele_1_count int(6) default NULL,
  allele_2_base char(1) default NULL,
  allele_2_count int(6) default NULL,
  count_N int(6) default NULL,
  used_coverage int(6) default NULL,
  based_on varchar(20) default NULL,
  repeat_type varchar(20) default NULL,
  exp_coverage double default NULL,
  coverage int(6) default NULL,
  PRIMARY KEY  (chromosome, position, ecotype),
  KEY chr_pos_index (chromosome,position),
  KEY ecotype_index (ecotype)
);



CREATE TABLE poly_cvp_unique (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  position int(9) unsigned NOT NULL default '0',
  reference char(1) default NULL,
  used_coverage int(6) default NULL,
  allele_1_base char(1) default NULL,
  allele_1_count int(6) default NULL,
  allele_2_base char(1) default NULL,
  allele_2_count int(6) default NULL,
  count_N int(6) default NULL,
  count_D int(6) default NULL,
  based_on varchar(20) default NULL,
  repeat_type varchar(20) default NULL,
  exp_coverage double default NULL,
  coverage int(6) default NULL,
  PRIMARY KEY  (chromosome, position, ecotype),
  KEY chr_pos_index (chromosome,position),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_cvp_duplication (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned NOT NULL default '0',
  length int(9) unsigned default NULL,
  rnv_count int(4) unsigned default NULL,
  obs_2_exp_cov double default NULL,
  PRIMARY KEY  (ecotype, chromosome, begin, end),
  KEY chr_pos_index (chromosome, begin, end),
  KEY ecotype_index (ecotype)
);

CREATE TABLE poly_cnv (
  ecotype varchar(20) NOT NULL default '',
  chromosome int(7) unsigned NOT NULL default '0',
  begin int(9) unsigned NOT NULL default '0',
  end int(9) unsigned default NULL,
  length int(9) unsigned default NULL,
  num_cvp_uniq int(4) unsigned default NULL,
  num_cvp_repetitive int(4) unsigned default NULL,
  num_repetitive int(9) unsigned default NULL,
  avg_ratio double default NULL
);
