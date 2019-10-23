
###############################################  mixed aligens with AOGC
#/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis
## /media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt
## /media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##Org
#analysis.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis"
#annotate.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Annotate"


##AN: start for novoaligned

analysis.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-11-11_AML/Analysis"
annotate.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-11-11_AML/Annotate"

## /media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chr1.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt
##/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-10-09_AML_CompleteGenomics_HC/Analysis/2014-10-09_AML_CompleteGenomics_HC.chrALL.Indel.ALL.ALL_GENOTYPES_analysis.txt

#################### use extension and names to select correct file in folder /Analysis
##Org
# project.extension<-"_analysis-maf-filtered.txt"
# project.name<-"2015-08-14_AML_mixedAligners." ## prefix for output file
#AN
project.extension<-"_analysis.txt"
project.name<-"2015-11-11_AML" ## prefix for output file
#change as per analysis: "type.project"
path.save<- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/case_control/"
type.project<-"case_control"

#fam<-c("HC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect


#fam<-c(".ALL.ALL_GENOTYPES") #  ALL or  c() ""-one project (the prefix of the summary files to collect ##Org
fam<- c("SUBSET2.ALL.ALL_GENOTYPES") ##AN
run.per.chromsome<-FALSE # TRue is doing per chromosome
use.genotype.recovery<-FALSE
force.recovery.model<-FALSE

#the.sample.sheet<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2013-02-27_AML_with_AOGCControl/BAM/TGCM-AML-combine_SampleSheet.csv"
#the.sample.sheet<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/CompleteGenomics/Chort_descriptions/Sequencing comparisons-ver 6.csv"
the.sample.sheet<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/sample_sheet.for.analysis.txt" ##org
#Use org, this file below does not have all the information required
#the.sample.sheet <- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/2015-11-11_AML.chr.ALL.Sample_Sheet_NEW.txt"  ##AN

related.or.bad.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/Analysis/related.or.bad.txt"
#pca.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/SNPs/pca_aml.txt" ##Org
pca.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/skat_test/pca_aml_novoaligned_march2016.txt.csv" ##AN

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)
remove.cols<-c()

#regions.file<-"/media/scratch2/AOGC-NGS/GFOS/gefos.seq/METHODS/0613-skatmeta-gefos/static/Homo_sapiens.GRCh37.70.protein_coding.genespace_boundaries.5k.split100k.txt"
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-FALSE ##

GATK.SB<-TRUE


a.label<-"CoVarRun.noControl.AML.regions"
maf.threshold.filter.to.use<-c(0.001,0.005,0.01,0.025)

sample.types<-c("Control","AML","AML-Child","AML-NotDiagnosis-Child","Asian-AML","Asian-AML-Child","Asian-AML-NotDiagnosis-Child","Asian-Control") ## these are the disease cohorts
case.control.classes<-c(0,1,9,9,9,9,9,9) ### Must be the same length as above

PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]



hw.target<-"Control"
high.missing.subset<-c("AML.filt","Control.filt")
cancer.group<-"AML"  ## used for no.genotypes.cancer etc
control.group<-"Control"
PD.group<-"PD" # if not defined in samples types this is not used


PD.group.types<-c("AML-Child","Asian-AML-Child","Asian-AML","AML-NotDiagnosis-Child","Asian-AML-NotDiagnosis-Child")
to.unwind.name<-"EVERYTHING"
snap.file.prefix<-"AML_novoalign_"

missing.targets<-c("AML.filt","Control.filt")
targets.analysis<-c("AML","Control") 


SnpStats.file.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats"
#SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"


#SnpStats.read.filter.file.suffix<-".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt$" ##Org
#SnpStats.stats.file.suffix<-".ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt$" ##org

SnpStats.read.filter.file.suffix<-"SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.read.filter.txt$" ##AN
SnpStats.stats.file.suffix<-"SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$" ##AN


##Org
# #SnpStats.bam.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
# SnpStats.bam.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt" ##*** May need updating
# SnpStats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
# SnpStats.stats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-08-14_AML_mixedAligners/snpstats/2015-08-14_AML_mixedAlignersHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"
# ###########################################################################

##AN
SnpStats.bam.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/2015-11-11-AML_snpstat/snpstat_input_novoalign_copy.txt" ##*** May need updating
##sample.rf
SnpStats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-11-11_AML/snpstats/2015-11-11_AMLnovoaligned.HC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.read.filter.txt"
##FAD
SnpStats.stats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-11-11_AML/snpstats/2015-11-11_AMLnovoaligned.HC.ALL.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt"
###########################################################################



a.indel.stat.use <- TRUE ## if false makes a dummy FALSE of have no snpstates
FAD.col<-"FAD" # usually "FAD" col in a.indel.stats used for alleleic depth AD
DUP.col<-"DUP" # usually "DUP"
#if(target.pheno.col %in% case.control){
covars<-c("PCA.1","PCA.2","PCA.3","PCA.4") #  covars<-c("1") # c("AGE_SCAN","PCA1","PCA2","PCA3","PCA4") #AGE_SCAN,PCA1,PCA2,PCA3,PCA4 #covars<-c("1")
#covars<-c("1")
cohorts.to.exclude.from.assoc<-c("SDDS","MODY","SKDP") #  table(pheno.ori[, "Project"])
subset.somatic<-FALSE # is true remove nomrmal gentypes
do.MAFS<-c(0.001) # do.MAFS<-c(0.001,0.01)
types<-c("coding","non-coding","sliding.window")

#####################################################

################################## Other input files needed - path required in file names
num.cores<-3 # use for filtering don't make too large else will exceed memory

gene.symbol.file.for.clusters<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Software/annovar/humandb/Gene_symbol_aliases.txt"
cluster.definition.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Final_FANC_clusters.csv" ## tab delimited with header
ensembleID.to.HGNC.map.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Genomes/hg19/ENSG_to_HGNC.txt" # tab delimited with header
coverage.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Mon_Jul_20_2015.txt" # tab delimited with header


## #/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt
## SnpStats.file.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Snp_read_filter"
## SnpStats.file.suffix<-"HC.2.BEST.chrALL.ACC_SUBSET.ALL.ALL_GENOTYPES_analysis.txt_snpstats.stats.txt$"
## #SnpStats.bam.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input.txt"
## SnpStats.bam.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/snpstats_input_aligner.txt"
## SnpStats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.read.filter.txt"
## SnpStats.stats.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/snp_stats/2015-03-16_AllAMLandLungHC.2.BEST.chrALL.ACC_SUBSET2.ALL.ALL_GENOTYPES_analysis-maf-filtered.txt_snpstats.stats.txt"

BWA.bad<-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/containamination loci.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)


###############
library(skatMeta)  ## ridge regression
#library(SKAT) ## skat method
library(GenomicFeatures)
library(HardyWeinberg)
library(Biostrings)

options(width=250,max.print=5000)

code.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
setwd(code.dir)

source("annotate_SNPs_subroutines.r")
source("Sliding_window_generation.r")
source("hwe.r")



contaminated.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/AOGC_exome_chip/Phenotypes/contaminated_AOGC_SEQ_samples.txt"
contaminated<-read.table(contaminated.file,header=F,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

#related.or.bad.file<- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/related.or.bad.txt"
related.or.bad<-read.table(related.or.bad.file,header=T,fill=TRUE,sep="\t",stringsAsFactors=FALSE)

if(SnpStats.bam.file!="" | exists(SnpStats.bam.file)) {
  #SnpStats.bam<-read.delim(SnpStats.bam.file,header=F,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE,check.names=FALSE) ##org
  SnpStats.bam<-read.delim(SnpStats.bam.file,header=F,sep=",",fill=TRUE,skip=0,stringsAsFactors=FALSE,check.names=FALSE) ##AN
  # SnpStats.bam[1:5,]
}




###########################
#/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Fri_Feb_06_2015.txt  latest
core.ann<-c("chr","start","end","REF","ALT","TYPE") # out put to annanlsys programs and need foe colun labels
dont.build.summary<-TRUE ##
GATK.SB<-TRUE
#max(BWA.bad[,"P.value"])




###################################### load ola aogc######################3
## load("/media/Bioinform-D/Research/annovar/humandb/aogc.count.data.RData")
## #print(colnames(indels))
## print(colnames(geno.aogc))
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70  ### this is where to add AOGC data INSERTS AFTER THE LOCATION
## ################################ add aogc

## ann<-readRDS("/media/scratch2/AOGC-NGS/Analysis/AOGC_sequnnce_LS/AOGC_sequence_10_LS_ANNOTATION.rds")




geneanno.DB<-c("refGene","knownGene","ensGene") # returns 2 extra columns
names(geneanno.DB)<-c("refGene","knownGene","ensGene")





##################################################### GEFOS - GENE LIST #####################################################
#mafs<-colnames(a.indel)[grepl("maf",colnames(a.indel))]
############################################# POPULATION MAF FILTER - PART A
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER

maf.threshold<-0.0  #  MAF threshold for annovar calling zero useful to get back all results !!do not modify!!

maf.threshold.filter.to.use<-sort(as.numeric(maf.threshold.filter.to.use))

filter.cols.novel.use<-c("PopFreqMax","NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHLBI_5400_ALL","NHLBI_5400_EUR","NHLBI_5400_AFR","1000genome","1000genome_asian","1000genome_mine","snp141","snp141_clinical","snp137","CG69","EUR_ASN_AFR_INDEL","AOGC-NGS_ALL","AOGC-NGS_ALL_OLD","Chinese") ##
#filter.cols.maf.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135","AOGC-NGS_ALL")
filter.cols.maf.use<-c("NHBLI_6500_ANNOVAR_ALL","NHBLI_6500_ALL","NHBLI_6500_EA","NHBLI_6500_AA","NHLBI_5400_ALL","1000genome","snp141","snp137","snp135")  # don't use AOGC if they are controls
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER



######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
############################################SET UP FNCTIONAL  FILTERS #######################################################

####################### MUTATION TYPE DEINITIONS 

possible.mutations<-c("frameshift substitution","nonframeshift substitution","downstream","frameshift deletion","frameshift insertion","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","ncRNA_splicing","ncRNA_UTR3","ncRNA_UTR5","ncRNA_UTR5;ncRNA_UTR3","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing","stopgain SNV","stoploss SNV","synonymous SNV","unknown","upstream","upstream;downstream","UTR3","UTR5","UTR5;UTR3","stopgain","stoploss")

interesting.coding.mutations<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","stopgain","stoploss")

interesting.mutations.use<-c("frameshift substitution","nonframeshift substitution","nonframeshift deletion","nonframeshift insertion","frameshift deletion","frameshift insertion","nonsynonymous SNV","stopgain SNV","stoploss SNV","splicing","ncRNA_exonic","stopgain","stoploss")

wanted.noncoding.subtypes<-c("miRNA","lincRNA") # filter by interesting to prefiler and vep.noncoding so dones get ncRNA intronic ::use gerp.score.threshold.low only these subtypes

interesting.to.prefilter<-c("UTR3","UTR5","UTR5;UTR3","snoRNA","snRNA","antisense","sense_intronic","ncRNA_exonic","ncRNA_splicing","intronic") #use gerp.score.threshold

extra.vep.annotations<-c("Uploaded_variation","Gene","Feature","Protein_position","Amino_acids")

vep.types<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")

## vep.types<-c( "not_assigned","stop_gained","stop_lost","stop_lost,NMD_transcript_variant","stop_gained,splice_region_variant,NMD_transcript_variant","initiator_codon_variant,splice_region_variant","splice_region_variant,3_prime_UTR_variant","stop_gained,NMD_transcript_variant","missense_variant,splice_region_variant","missense_variant","splice_acceptor_variant","splice_acceptor_variant,nc_transcript_variant","splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant","splice_donor_variant,nc_transcript_variant","splice_region_variant,intron_variant,NMD_transcript_variant","splice_donor_variant","splice_region_variant","splice_region_variant,5_prime_UTR_variant","splice_region_variant,synonymous_variant","splice_region_variant,intron_variant,nc_transcript_variant","splice_region_variant,non_coding_exon_variant,nc_transcript_variant","missense_variant,NMD_transcript_variant","splice_region_variant,intron_variant","NMD_transcript_variant","intron_variant,NMD_transcript_variant","mature_miRNA_variant","5_prime_UTR_variant","5_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant,nc_transcript_variant","3_prime_UTR_variant,NMD_transcript_variant","non_coding_exon_variant","TF_binding_site_variant","intron_variant,nc_transcript_variant","synonymous_variant,NMD_transcript_variant","3_prime_UTR_variant","regulatory_region_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","intron_variant","synonymous_variant")


vep.coding<-c("not_assigned","stop_gained","stop_lost","missense_variant","splice_acceptor_variant","splice_donor_variant","splice_region_variant","initiator_codon_variant","stop_retained_variant","incomplete_terminal_codon_variant","frameshift_variant","inframe_deletion","inframe_insertion")

vep.noncoding<-c("5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","NC_stop_gained","NC_stop_lost","NC_splice_acceptor_variant","NC_splice_donor_variant","NC_splice_region_variant","NC_initiator_codon_variant","NC_stop_retained_variant","NC_non_coding_exon_variant","NC_incomplete_terminal_codon_variant","NC_3_prime_UTR_variant","mature_miRNA_variant","NC_5_prime_UTR_variant","TF_binding_site_variant","TFBS_ablation","TFBS_amplification","regulatory_region_variant","intron_variant")              

vep.unwanted<-c("intron_variant","NC_intron_variant","synonymous_variant","coding_sequence_variant","NC_synonymous_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","NC_intergenic_variant","NMD_transcript_variant","nc_transcript_variant","NC_nc_transcript_variant","feature_truncation","feature_elongation")



missense.variant<-c("nonsynonymous SNV","missense_variant")
synonymous.variant<-c("synonymous SNV","synonymous_variant")

hwe.control.threshold<-1e-8
gerp.score.threshold.high<-2.5 # gerp score >= will be included
gerp.score.threshold.low<-2.0 # gerp score >= will be included
gerp.score.threshold.unknown<-0


#generic.filter.DB

interesting.mutations<-interesting.mutations.use


if(GATK.SB){
  global.quality.labs<-c("QUAL","QD","HRun","SB","FILTER","FILTER","PolyPhen.scores","PolyPhen.scores","SIFT.scores","mut.taster::score","phylo::score","PolyPhen.desc","SIFT.desc","GERP::score","GERP::score","GERP::score","MAF.ALL","MAF.HIGH","MAF.LOW","TYPE") ### THESE ARE THE COLUMN LABELS IN THE DATA these become the "good.qual" filter 
  global.quality.names<-c("QUAL","QD","HRun","SB","FILTER_PASS","FILTER_100","PolyPhen.low","PolyPhen.high","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","GERP.low","GERP.unknown","MAF.ALL","MAF.HIGH","MAF.LOW","flat") ### THESE ARE THE COLUMN LABELS IN the quality.filter TABLE
  #global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",2,2,0.25,0.25,0.25)
  global.quality.cut<-c(50,0.5,5,1,"PASS","TruthSensitivityTranche99.90to100.00",0.1,0.4,0.4,0.4,0.4,"damaging","deleterious",gerp.score.threshold.high,gerp.score.threshold.low,gerp.score.threshold.unknown,0.25,0.25,0.25,"flat")
  global.quality.type<-c("numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","factor","factor","numeric","numeric","numeric","numeric","numeric","numeric","factor")
  global.quality.dirn<-c("greater","greater","less","less","exact","exact","greater","greater","greater","greater","greater","ends_with","exact","greater","greater","exact","greater","greater","greater","ends_with")
}else{
  global.quality.labs<-c("QUAL","QD","HRun","FS","FILTER","FILTER") ### these become the "good.qual" filter
  global.quality.names<-c("QUAL","QD","HRun","FS","FILTER_PASS","FILTER_100")
  global.quality.cut<-c(50,0.5,5,60,"PASS","TruthSensitivityTranche99.90to100.00")
  global.quality.type<-c("numeric","numeric","numeric","numeric","factor","factor")
  global.quality.dirn<-c("greater","greater","less","less","exact","exact")
}
names(global.quality.cut)<-global.quality.labs
names(global.quality.dirn)<-global.quality.labs
names(global.quality.type)<-global.quality.names
global.labs<-unique(global.quality.labs)

global.quality.cut
global.quality.dirn
global.quality.type

quality.cut<-global.quality.cut
quality.type<-global.quality.type
quality.dirn<-global.quality.dirn
######################################################################################################################

#"/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Wed_Feb_04_2015.txt"  # leo coverage
seq.type.file<-coverage.file  ## Troels-covergae prior to top up
seq.type<-read.delim(seq.type.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
coverage<-seq.type # read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/QC for all samples summary/Coverage_QC/QC_stat_SAMPLE_Tue_Oct_14_2014.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)

#to.fix<-(coverage[,"Project"]=="SDDS")
#coverage[to.fix,"Sample"]<-gsub("^0+","", coverage[to.fix,"Sample"])  ## SDDS project has proceeding zeros


seq.type[1:5,]
seq.type.sub<-seq.type[seq.type[,"Project"]=="AOGC-NGS",]
nim.samples<-seq.type.sub[seq.type.sub[,"Capture.Method"]=="TruD:NimX","Sample"]
ill.samples<-seq.type.sub[seq.type.sub[,"Capture.Method"]=="TruD:TruX","Sample"]

nim.samples<-paste(nim.samples,"GT",sep=".")
ill.samples<-paste(ill.samples,"GT",sep=".")


length(nim.samples)
length(ill.samples)

############################################ GET SAMPLE INFO and PHENOTYPES ####################################################

## ########check

print(paste0("The samples sheet is: ",the.sample.sheet))

sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) #,check.names=FALSE )
sample.sheet.full[1:5,1:5]
colnames(sample.sheet.full)
dim(sample.sheet.full)
table(sample.sheet.full[,"SampleProject"])
table(sample.sheet.full[,"Sex"])
##### fix 0 and 9 for missing to NA


########################################## ADD PCA data


if(exists("pca.file")){
  #pca<-read.delim(pca.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE) ##org
  pca<-read.delim(pca.file,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE) ##AN
  
  pca[1:5,]
  sample.sheet.full[1:5,]
  
  posns<-match(sample.sheet.full[,"ParticipantCode"],pca[,"Sample"])
  missing<-is.na(posns)
  sum(missing)
  print(paste0("samples missing from pca file that in the the sample.sheet ",sample.sheet.full[missing,"ParticipantCode"]))
  
  extra<-pca[posns,]
  sample.sheet.full<-cbind(sample.sheet.full,extra)
  print(sample.sheet.full[1:2,])
}

################################################# coverage check
################################################# coverage check
######### Check and fix the sample sheet


dim(coverage)
coverage[1:5,]
unique(coverage[,"Project"])
aml.projects<-c("AML-exome","RSGB_AML","AMAS", "TGCM-AML")
aml.samples<-coverage[coverage[,"Project"] %in% aml.projects,"Sample"]
length(aml.samples)
aml.have<-sample.sheet.full[,"ParticipantCode"] %in% aml.samples
table(sample.sheet.full[,"SampleProject"])
table(sample.sheet.full[aml.have,"SampleProject"])
table(sample.sheet.full[!aml.have,"SampleProject"])

table(coverage[, "Capture.Method"])

#SnpStats.bam[1:5,]

if(exists("coverage")){
  
  posns<-match(sample.sheet.full[,"ParticipantCode"], coverage[,"Sample"])
  missing<-is.na(posns)
  sum(missing)
  seq.project<-coverage[posns, c("Project","Capture.Method")]
  capture<-coverage[posns, "Capture.Method"]
  table(capture)
  capture<-strsplit(capture,split=";")
  lengths<-unlist(lapply(capture,length))
  capture[lengths==0]<-"unknown"
  
  capture<-unlist(lapply(capture,function(x) { paste(unique(x),collapse=";")})) #x[(is.na(x) | x=="")]<-"NONE";
  table(capture)
  #cbind(capture,seq.project)
  
  sample.sheet.full<-cbind(sample.sheet.full,capture,seq.project)
}
table(sample.sheet.full[,"capture"])
table(sample.sheet.full[,"capture"])


if(exists("SnpStats.bam")){
  posns<-match(sample.sheet.full[,"ParticipantCode"],SnpStats.bam[,2])
  missing<-is.na(posns)
  sum(missing)
  Aligner<-SnpStats.bam[posns,3]
  sample.sheet.full[is.na(Aligner),"ParticipantCode"] ### these wree all aogc
  Aligner[is.na(Aligner)]<-"novoalign"
  sample.sheet.full<-cbind(sample.sheet.full,Aligner)
}
table(sample.sheet.full[,"Aligner"])
table(sample.sheet.full[,"Aligner"])
colnames(sample.sheet.full)
## bwa     novoalign novoalign;bwa 
##  81           557            11
################################################# Build snpstats file
################################################# coverage snpstats file





## a.control<-grepl("Control",sample.sheet.full[,"SampleProject"])
## a.case<-grepl("AML",sample.sheet.full[,"SampleProject"])
## AFF_STATUS<-rep(NA,times=dim(sample.sheet.full)[1])
## sum(a.case & a.control) # should be 0
## AFF_STATUS[a.case]<-2
## AFF_STATUS[a.control]<-1

## #sample.sheet.full[aml.have & a.control,]

## posns<-match(sample.sheet.full[,"ParticipantCode"], coverage[,"Sample"])
## missing<-is.na(posns)
## sum(missing)

## sample.sheet.full.coverage<-coverage[posns,]
## sample.sheet.full.coverage[1:5,]

## SnpStats.input.file<-cbind(paste("/mnt/UQCCG/Sequencing/Projects/",sample.sheet.full.coverage[,"Project"],"/BAM/",sample.sheet.full.coverage[,"BAM"],sep=""),sample.sheet.full.coverage[,"ID"],AFF_STATUS,sample.sheet.full[,"ParticipantCode"])
## SnpStats.input.file[1:5,]

## colnames(SnpStats.input.file)<-c("BAM","SAMPLE_ID","AFF_STATUS","ParticipantCode")

## all.bams<-SnpStats.input.file[,"SAMPLE_ID"]


## mult.bams<-strsplit(all.bams,split=";")
## num.mult.bams<-unlist(lapply(mult.bams,length))
## first.bam<-unlist(lapply(mult.bams,function(x) x[1]))
## #has.mult.bams
## SnpStats.input.file[,"SAMPLE_ID"]<-first.bam
## SnpStats.input.file<-cbind(SnpStats.input.file,num.mult.bams,all.bams)
## SnpStats.input.file[1:5,]


## getwd()
## write.table(SnpStats.input.file,file="SnpStats.input.file.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## ## /media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/2015-03-16_AllAMLandLung.chr13.ALL.ALL_GENOTYPES_analysis.txt



##  [1] "860"  "861"  "862"  "2529" "2530" "2531" "2567" "2688" "2689" "2851" "2924" "2944" "3198" "3199" "3200" "3204" "3255" "3256" "3257" "3258" "3259" "3342" "3343" "3344" "3345" "2568" "2569" "2669" "2670" "2671" "3040" "3041" "3169" "3170" "3211"
## [36] "3212" "3282" "3287" "3288" "3289" "3291" "3293" "3294" "3295" "3309"


################################################# coverage check
################################################# coverage check



########################################################### define bad samples
#bad.samples<-sample.sheet.full[ (sample.sheet.full[,"SampleProject"]!="AML" & sample.sheet.full[,"SampleProject"]!="Control"  )   ,"ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"
#table(sample.sheet.full[ !(sample.sheet.full[,"ParticipantCode"] %in% bad.samples) ,"SampleProject"])

## AML Control 
## 136     363
## tail(sample.sheet.full)
## test<-((sample.sheet.full[,"Race"]!="W" | is.na(sample.sheet.full[,"Race"])) & sample.sheet.full[,"AffectionStatus"]==2) 
## sample.sheet.full[test,]
## bad.samples<-sample.sheet.full[test,"ParticipantCode"] # "1"  "9"  "63" "74" "83" "84" "99"

contaminated.aogc<-sample.sheet.full[sample.sheet.full[,"ParticipantCode"] %in% contaminated[,1],"ParticipantCode"]
related.or.bad[,1]
bad.samples<-unique(c(contaminated.aogc,related.or.bad[,1]))

table(sample.sheet.full[ !(sample.sheet.full[,"ParticipantCode"] %in% bad.samples) ,"SampleProject"])
################################################################################









pheno.types<-c("SampleProject") ## vales is column header
names(pheno.types)<-c("SampleProject") ### name is output columns

## case.control.classes<-c(0,1,0)
## names(case.control.classes)<-c("Control","NMD","AOGC")
## case.control.classes


case.control<-c("SampleProject")
sample.types
# case.control.classes<-c(0,1,9,9,9,9,9,9)
names(case.control.classes)<-sample.types
case.control.classes
# ib<-1
for(ib in 1:length(case.control)){
  if(!(case.control[ib] %in% colnames(sample.sheet.full))){next}
  ##Ask paul!!!
  sample.sheet.full[(  !(sample.sheet.full[,case.control[ib]] %in% names(case.control.classes))  |  is.na(sample.sheet.full[,case.control[ib]]) | sample.sheet.full[,case.control[ib]]==0 | sample.sheet.full[,case.control[ib]]==9)  ,case.control[ib]]<-NA
}

table(sample.sheet.full$SampleProject)

## AML                    AML-Child       AML-NotDiagnosis-Child                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control                      Control 
## 136                           22                           17                           11                            2                            4                           23                          363


#tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)
colnames(sample.sheet.full)[colnames(sample.sheet.full)=="ParticipantCode"]<-"SAMPLE"
control.samples<-{}
all.samples<-sample.sheet.full[,"SAMPLE"]

length(all.samples)

## o.remove.all<-expand.labels.to.samples(remove.from.all.samples,all.samples)
## to.remove.samples<-unique(to.remove.all)
## remove.cols<-unique(c(remove.cols,to.remove.samples))







#### test fam list
files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
toString( unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   ))) )
toString( files)
files
####


#############################################################################################################

######################################### Predefined variables required
##################################################################################


#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)
getwd()

files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
if(fam=="ALL" | fam=="All" | fam=="all"){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam
## ifam<-2

## geno.aogc<-read.delim("/media/scratch2/AOGC-NGS/Analysis/AOGC-Genotyping.output.AOGC_ALL.geno.all.txt",header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
## setwd("/media/Bioinform-D/Research/annovar/humandb/")
## save(list=c("geno.aogc","use.key"),file="aogc.count.data.RData")


#load("/media/Bioinform-D/Research/annovar/humandb/aogc.count.data.RData")
#print(colnames(indels))
## print(colnames(geno.aogc))
## use.key<-build.key(geno.aogc,core.ann)
## insert.location<-70 
### this is where to add AOGC data
################################ add aogc

## project.files
######################################################################################################
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##############    ifam begins    ####################################################################
######################################################################################################
##AN: commented out : for(ifam in 1:length(fam)){
ifam<-1
#for(ifam in 1:length(fam)){
  
  the.extension<-paste(fam[ifam],project.extension,"$",sep="")
  project.files<-files[grepl(the.extension ,files)]
  print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]
  
  indels<-{}
  the.col<-{}
  project.files
  # ichr<-13
  project.files
  
  
  if(run.per.chromsome){
    project.files<-project.files[!grepl("chrALL",project.files)]
    
    if(length(project.files)!=24){
      print("########################################### WARNING #################################")
      print("less that 24 chromosomes detected")
      print(fam[ifam])
      print("########################################### WARNING #################################") 
    }
  }
  project.files
  
  
  ######################################################################################################
  ##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ##############    ichr begins    ####################################################################
  ######################################################################################################
  ##AN: commented out :for(ichr in 1:length(project.files)){
  
  ichr<-1
  
 # for(ichr in 1:length(project.files)){
    


print("START")
setwd(analysis.dir)
################## fast read ###########  load("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Analysis/LEO_Pharma_wanted_with_synon_ALL_0.01_muts_Feb23")
column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
print(project.files[ichr])
the.chr<-gsub(the.extension,"",project.files[ichr])
the.chr<-gsub(project.name,"",the.chr)
print(the.chr)
## load("2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt_SUBSET.RData")
## column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
## length(column.labels)
## dim(a.indel)
## colnames(a.indel)<-column.labels
## load("2015-03-16_AllAMLandLung.BEST.chrALL.ACC_0.025.ALL.ALL_GENOTYPES_analysis.txt.RData")
column.labels<-colnames(a.indel)
a.indel<-as.matrix(a.indel)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
## anted.genes<-c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","WT1","MLL","KRAS","FLT3","IDH2","IDH1","TP53","KIT","NPM1","JAK2","DNMT3A","TET2","RUNX1","NRAS","CEBPA","PTPN11","U2AF1","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","HNRNPK","FANCA","FANCB","FANCC","FANCD1","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","C1orf86","C19orf40","C17orf70","SLX1","MUS81","ERCC1","FAN1","EME1","EME2","MRE11A","NBN1","RAD50","FAND1","BRCA1","BARD1","RAD51","RAD51B","RAD51D","XRCC2","XRCC3","RMI1","RMI2","BLM","TOP3A","RPA1","RPA2","RPA3","ATM","ATR","ATRIP","CHECK1","RAD9A","RAD17","CS","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","ACLY","ACO1","OGDH","ACO2","PC","PCK1","PCK2","PDHA1","PDHA2","PDHB","OGDHL","SDHA","SDHB","SDHC","SDHD","SUCLG2","SUCLG1","SUCLA2")
## wanted.genes2<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/Gene_Lists/mam_new/ALL_genes.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
## wanted.genes<-unique(c(wanted.genes,wanted.genes2[,1]))


## ALDH1B1:c.1132C>T (p.Gln378Ter) chr9

## ALDH2:c.1510G>A (p.Glu504Lys) chr12

###########################################################################test a single gene here
###########################################################################test a single gene here
###########################################################################test a single gene here
###########################################################################test a single gene here
###########################################################################test a single gene here
###########################################################################test a single gene here 
## wanted.genes<-c("RHEBL1")
## wanted.genes<-c("CPB1")
## wanted.genes<-c("FRY")
## wanted.genes<-c("WDR1")

## out.name<-toString(wanted.genes)

## ## wanted<-grepl("ALDH1B1",a.indel[,"Gene.Names"])
## ## aa.wanted <-grepl("1132",a.indel[,"refGene::type"])

## ## wanted<-grepl("ALDH2",a.indel[,"Gene.Names"])
## ## aa.wanted <-grepl("1510",a.indel[,"refGene::type"])

##  wanted<- a.indel[,"Gene.Names"] %in% wanted.genes

## a.indel[wanted,"refGene::type"]
## a.indel[wanted,c("Gene.Names","ensGene::type","ID::maf")]



## sum(wanted)
## sum(aa.wanted & wanted)
## a.indel[aa.wanted & wanted,1:20]

## a.indel<-a.indel[wanted,]
## aa.wanted <-grepl("1132",a.indel[,"refGene::type"])
## a.indel[aa.wanted,c("Gene.Names","ensGene::type","ID::maf")]

## key[aa.wanted]

## #### get the summary.geno.extra at line 1660
## summary.geno.extra[,grepl("GENO",colnames(summary.geno.extra))]
## summary.geno.extra[aa.wanted,grepl("GENO",colnames(summary.geno.extra))]
## colnames(a.indel)

## pheno<-pheno.ori

## the.samples.use
## ###################################
## #PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]
## genotypes.ex<-a.indel[, grepl(".GT$",colnames(a.indel)) ]
##  dim(genotypes.ex)
## genotypes.ex<-t(genotypes.ex)

## genotypes.ex[genotypes.ex=="NA"]<-NA
## genotypes.ex[genotypes.ex=="0/0"]<-0
## genotypes.ex[genotypes.ex=="0/1"]<-1
## genotypes.ex[genotypes.ex=="1/1"]<-2

## rownames(genotypes.ex)<-gsub(".GT","",rownames(genotypes.ex))

## dim(genotypes.ex)
## dim(genotypes.ex)
## ###############################################
## posns<-match(pheno.ori[,"SAMPLE"],rownames(genotypes.ex))
## missing<-is.na(posns)
## rownames(genotypes.ex)[posns[missing]]
## genotypes.ex<-genotypes.ex[posns[!missing],]
## dim(genotypes.ex)
## dim(pheno.ori)
## sum(rownames(genotypes.ex)!=pheno.ori[,"SAMPLE"])
## sum(pheno.ori[,"SAMPLE"]!=rownames(genotypes.ex))

## muts.in.PD<-apply(genotypes.ex[pheno.ori[,"PD"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.cases<-apply(genotypes.ex[pheno.ori[,"AML"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
## muts.in.controls<-apply(genotypes.ex[pheno.ori[,"Control"],],2,function(x) { paste(names(x)[x!=0 & !is.na(x)],collapse=",")})
##  ann.cols<-c("chr","start","end","REF","ALT","TYPE","refGene::type","knownGene::type","ensGene::type","Gene.Names","Genes.mentioned.at.ASH","refGene::location","knownGene::location","ensGene::location","OMIM (Gene::Status::OMIM::description::disease)","Consequence.Embl","Uploaded_variation.Embl","Gene.Embl","Feature.Embl", "Protein_position.Embl", "Amino_acids.Embl" , "ensGene::type","ID::maf","FILTER")

## output<-cbind(a.indel[,ann.cols],muts.in.cases,muts.in.PD,muts.in.controls,summary.geno.extra[,colnames(summary.geno.extra)[grep("^GENO",colnames(summary.geno.extra))]])
## output[aa.wanted,]
## getwd()
## write.table(output,file=paste0(out.name,"_genotype.summary.txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## #write.table(a.indel,file="ALDH2_genotypes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


## ## write.table(output,file="ALDH2_genotype.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## ## write.table(a.indel,file="ALDH2_genotypes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

## ## write.table(output,file="ALDH1B1_genotype.summary.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
## ## write.table(a.indel,file="ALDH1B1_genotypes.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


## output.1<-output
###########################################################################END test a single gene here
###########################################################################END test a single gene here
###########################################################################END test a single gene here
###########################################################################END test a single gene here

##########################################################################
########################### START ADDITIONAL FILTERING " : bad.qual.locations############################
##########################################################################
##########################################################################
## extra filtering
#SNP_STATS:version1
#FILTER_SUMMARY:alt_mismatch5_proportion_gt_0.8;Strand_bias_proportion_gt_0.9;Read_end3_proportion_gt_0.8;High_not_called_alts_gte_3
#SUMMARY_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;prop_mismatch;prop_fstrand;prop_rstrand;prop_read_end
#SUMMARY_NOT_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;high_alt_count
#SAMPLE.RF:GATK_called_snp;ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count

#####################################read filter file

if(SnpStats.file.dir!="" | SnpStats.file!=""){
if(run.per.chromsome){
    
files<-dir(SnpStats.file.dir)
files<-files[grepl(paste(the.chr,SnpStats.read.filter.file.suffix,sep=""),files)]
files
i<-1
setwd(SnpStats.file.dir)
for (i in 1:length(files)){
a.filt<-read.delim(files[i],header=T,sep="\t",fill=TRUE,skip=5,stringsAsFactors=FALSE)
if(i==1){
  filt<-a.filt
}else{
  filt<-rbind(filt,a.filt)
}

}

}else{ ## one names file not via chromosomes
  print(paste0("Reading SNPstats file",SnpStats.file))
  filt<-read.table(SnpStats.file,header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE)
  filt<-filt[,1:11]
} ## per chromosome

filt.key<-build.key(filt,core.ann,add.chr.label=TRUE)
setwd(analysis.dir)
# colnames(filt)
pass.filt<-strsplit(filt[,"FILTER_SUMMARY"],split=";")
pass.filt<-unlist(lapply(pass.filt,function(x) {sum(as.logical(x[1:3]))==0}))
#pass.filt[173713]
# sum(!pass.filt)
#filt.key[173713]
snp.fail.filt<-filt.key[!pass.filt]
# snp.fail.filt[1:5]

## AN:if - else

}else{ ### no snpstats file
    snp.fail.filt<-rep(TRUE,times=dim(a.indel)[1])
    extra<-matrix(data=NA,nrow=dim(a.indel)[1],ncol=length(c("FILTER_SUMMARY", "SUMMARY_CALLED", "SUMMARY_NOT_CALLED")))
    filt<-cbind(a.indel[,c(core.ann)],extra)
    colnames(filt)<-c(core.ann,c("FILTER_SUMMARY", "SUMMARY_CALLED", "SUMMARY_NOT_CALLED"))
    filt.key<-build.key(filt,core.ann,add.chr.label=FALSE)
}  # (SnpStats.file.dir!="" | SnpStats.file!="")


##################################### stats  file
if(SnpStats.file.dir!="" | SnpStats.file!=""){
    
if(run.per.chromsome){
    
files<-dir(SnpStats.file.dir)
files<-files[grepl(paste(the.chr,SnpStats.stats.file.suffix,sep=""),files)]
files
i<-1
setwd(SnpStats.file.dir)
for (i in 1:length(files)){
one.indel.stats<-read.table(files[i],header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
if(i==1){
  a.indel.stats<-one.indel.stats
}else{
  a.indel.stats<-rbind(a.indel.stats,one.indel.stats)
}

}

}else{ ## one names file not via chromosomes
    a.indel.stats<-read.table(SnpStats.stats.file,header=T,skip=5,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
}

}else{ ## no a.indel.stats so build a false one
    
    a.indel.stats<-a.indel[,c(core.ann,paste(all.possible.samples,FAD.col,sep="."))]
}# (SnpStats.file.dir!="" | SnpStats.file!="")




if(exists("one.indel.stats")){rm(one.indel.stats)}

   
if(exists("a.indel.stats")){
# a.indel.stats[1:5,1:12]
dim(a.indel.stats)
key.stats<-build.key(a.indel.stats,core.ann,add.chr.label=TRUE)
key.stats[1:5]
}

   
   
setwd(analysis.dir)   
## all.possible.samples.stats<-gsub(".FAD$","",colnames(a.indel.stats)[grep(".FAD$",colnames(a.indel.stats))],perl=TRUE)
## all.possible.samples.stats
## length(all.possible.samples)
## length(all.possible.samples.stats)
## all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]
##########################################################################
########################### END ADDITIONAL FILTERING ############################
##########################################################################

   

all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
all.possible.samples.stats<-gsub(".FAD$","",colnames(a.indel.stats)[grep(".FAD$",colnames(a.indel.stats))],perl=TRUE)
length(all.possible.samples)
length(all.possible.samples.stats)

print(paste0("Samples missing in genotyping vs Snpstats:",all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]))
#print(all.possible.samples[!(all.possible.samples %in% all.possible.samples.stats)]) #

pheno.types



the.chr<-a.indel[1,"chr"]
print(paste("Doing Chromosome ",the.chr))
print(table(a.indel[,"chr"]))

if(!grepl("^chr",the.chr)){
a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}
key<-build.key(a.indel,core.ann)
rownames(a.indel)<-key
length(key)
   length(key.stats)
key[1:5]
key.stats[1:5]
####### REMOVE BAD SAMPLES
## key.stats<-build.key(a.indel.stats,core.ann[1:5],add.chr=TRUE)
## key<-build.key(a.indel,core.ann[1:5])
## length(key)
## length(key.stats)
## sum(key!=key.stats)
## cbind(key,key.stats)[key!=key.stats,][1:10,]


posns<-match(key,key.stats)
missing<-is.na(posns)
sum(missing)

## an.indel<-grepl("^indel",a.indel[,"TYPE"])
## a.indel[missing & !an.indel ,core.ann][1:10,]


a.indel.stats<-a.indel.stats[posns,]

rownames(a.indel.stats)<-key

dim(a.indel)
dim(a.indel.stats)


   
############################################## remove any samples that fail QC


if(length(bad.samples)>0){

bad.samples.labels<-expand.labels.to.samples(bad.samples,c("GT","AD","DP","GQ"),paste.after=TRUE)
if(length(bad.samples.labels)>1){
a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% bad.samples.labels)]]
}

bad.samples.labels<-expand.labels.to.samples(bad.samples,c("FAD","TAD","DUP"),paste.after=TRUE)
if(length(bad.samples.labels)>1){
a.indel.stats<-a.indel.stats[,colnames(a.indel.stats)[!(colnames(a.indel.stats) %in% bad.samples.labels)]]
}

} # have bad samples so remove


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


########################################################################################################################
####################################### EXTRA snp filtering ##################################################
########################alignemnt issues files
####################### these are position from the leo pharma project

posns<-match(key,BWA.bad[,"gene"])
poss.model<-BWA.bad[posns, c("gene","P.value")]
pass.possonian.control.model<- poss.model[,"P.value"]> 2*(pnorm(abs(6), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))  | is.na( poss.model[,"P.value"])

length(pass.possonian.control.model)  ## these are known BWA problem sites
#length(pass)
sum(!pass.possonian.control.model)




########################alignement issues files from filt array

bad.qual.locations<-key %in% snp.fail.filt
sum(bad.qual.locations)
length(bad.qual.locations)
#length(pass)

####################################################################### redo the populatiom fileter excluding AOGC
############################################# POPULATION MAF FILTER 
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER
############################################# POPULATION MAF FILTER


maf.threshold.filter<-maf.threshold.filter.to.use
# > maf.threshold.filter
# [1] 0.001 0.005 0.010 0.025

filter.cols.maf<-filter.cols.maf.use
# > filter.cols.maf
# [1] "NHBLI_6500_ANNOVAR_ALL" "NHBLI_6500_ALL"         "NHBLI_6500_EA"          "NHBLI_6500_AA"          "NHLBI_5400_ALL"         "1000genome"             "snp141"                 "snp137"                 "snp135"   

filter.cols.novel<-filter.cols.novel.use
# [1] "PopFreqMax"             "NHBLI_6500_ANNOVAR_ALL" "NHBLI_6500_ALL"         "NHLBI_5400_ALL"         "NHLBI_5400_EUR"         "NHLBI_5400_AFR"         "1000genome"             "1000genome_asian"       "1000genome_mine"       
# [10] "snp141"                 "snp141_clinical"        "snp137"                 "CG69"                   "EUR_ASN_AFR_INDEL"      "AOGC-NGS_ALL"           "AOGC-NGS_ALL_OLD"       "Chinese"     
#############################get the filters so have minor allele frequencies:

the.combined.DBs<-gsub("::maf$","",colnames(a.indel)[grep("::maf$",colnames(a.indel))]) ## get maf databases
the.combined.DBs<-the.combined.DBs[!(the.combined.DBs %in% c("ID","snp132"))] ## remove artifacts
names.the.combined.DBs<-the.combined.DBs
#the.combined.DBs<-paste(the.combined.DBs,"::maf",sep="")
names(the.combined.DBs)<-names.the.combined.DBs

 present.cols<-filter.cols.maf %in% names(the.combined.DBs)
if(sum(!present.cols)>0){  ### asked to filter on columns that don't exist
  print("WARNING asked to filter using MAF with databaes that have not been annotated")
  print(paste("missing MAF Db:",toString(filter.cols.maf[!present.cols]),"-> ignored",sep=" "))
  filter.cols.maf<-filter.cols.maf[present.cols]
}

 present.cols<-filter.cols.novel %in% names(the.combined.DBs)
if(sum(!present.cols)>0){  ### asked to filter on columns that don't exist
  print("WARNING asked to filter using NOVEL with databaes that have not been annotated")
  print(paste("missing MAF Db:",toString(filter.cols.novel[!present.cols]),"-> ignored",sep=" "))
  filter.cols.novel<-filter.cols.novel[present.cols]
}

print(filter.cols.maf)
print(filter.cols.novel)

maf.target<-the.combined.DBs[names(the.combined.DBs) %in% filter.cols.maf] # the.combined.DBs
maf.target

for(i in 1:length(maf.threshold.filter)){
a.filter<-paste(filter.cols.maf,"::maf-filter::",maf.threshold.filter[i],sep="")
print(a.filter)
 assign(paste("filter.cols.maf",maf.threshold.filter[i],sep="."),value=a.filter)
}

maf.threshold.filter.all<-c(0.0,maf.threshold.filter)
assign(paste("filter.cols.maf",maf.threshold.filter.all[1],sep="."),value=paste(filter.cols.novel,"::maf",sep=""))
print(eval(as.name(paste("filter.cols.maf",maf.threshold.filter.all[1],sep=".")) ) )

   
#  a.indel[1:5,paste(the.combined.DBs,"::maf",sep="")]


   
################################################################################AN:Population based filtering based on minor allele frequency maf.threshold.filter where new columns are added in target.table if the maf is lesser than maf.threshold.filter[k]
## THE FILTER TABLE IS MADE FROM : the.combined.DBs<-c(filter.DB,generic.filter.DB,function.filter.DB)  only thes DBs can be used for MAF filtereing
target.table<-a.indel[,paste(the.combined.DBs,"::maf",sep="")]
#target.table[1:5,]
extra<-matrix(data=NA,nrow=dim(target.table)[1],ncol=length(maf.target)*length(maf.threshold.filter))

colnames(extra)<-expand.labels.to.samples.complex(maf.target,maf.threshold.filter,paste.after=TRUE,seperator="::maf-filter::")
#extra[1:5,]

target.table<-cbind(target.table,extra)
rm(extra)
####### generate maf.target appreded to the allel frequency
maf.threshold.filter  # filters to apply
maf.target  #<-the.combined.DBs[names(the.combined.DBs) %in% all.filter.cols.maf] # the.combined.DBs

k<-2
i<-2


for(k in 1:length(maf.threshold.filter)){ #maf.threshold.filter.all contains novel "0"
for(i in 1:length(maf.target)){
  an.inversion<-(target.table[, paste(maf.target[i],"maf",sep="::") ] > (1-maf.threshold.filter[k]) )  # 0.999 on a 0.01 filter is also a hit
  an.inversion[is.na(an.inversion)]<-FALSE
#  a.test<-(target.table[,paste(maf.target[i],"maf",sep="::")] < maf.threshold.filter[k]) # normal test for minot allele
  target.table[,paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]<-(target.table[,paste(maf.target[i],"maf",sep="::")] < maf.threshold.filter[k]) | an.inversion  ## incase ref and alt are inverted
  target.table[is.na(target.table[,paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]),paste(maf.target[i],"maf-filter",maf.threshold.filter[k],sep="::")]<-TRUE # if not found in database then assume has small allele frequency
#  ref.inversion[,k ]<-ref.inversion[,k] &   an.inversion
}
print(paste("Done ",maf.threshold.filter[k],sep=""))
    }

#target.table[1:5,]
###################


################################# Get allele frequency table
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies
maf.lt.all<-data.frame(key=key,stringsAsFactors=FALSE)

  
  imaf<-1
  for(imaf in 1:length(maf.threshold.filter.all)){
  a.filter.cols.maf<-eval(as.name( paste("filter.cols.maf",maf.threshold.filter.all[imaf],sep=".")  ))
  maf.lt<-rep(TRUE,times=dim(target.table)[1])

  imaff<-1
  for(imaff in 1:length(a.filter.cols.maf)){
    if(maf.threshold.filter.all[imaf]==0){ # different case for novel NOT found (rather than less than)
  #     maf.lt<-maf.lt & !target.table[,a.filter.cols.maf[imaff]]  ### this is correct as using the variable string name derived from filter.cols.novel
        maf.lt<-maf.lt & ( is.na(target.table[,a.filter.cols.maf[imaff]])  | target.table[,a.filter.cols.maf[imaff]]=="NA" )
     }else{
    maf.lt<-maf.lt & as.logical(target.table[,a.filter.cols.maf[imaff]])
  }
  }
# filtered<-maf.lt & wanted.muts.fil
 maf.lt.all<-cbind(maf.lt.all,maf.lt)
}

   
if(dim(maf.lt.all)[2]>1){ colnames(maf.lt.all)<-c("key",paste("MAF.lt:",maf.threshold.filter.all,sep=""))}
maf.lt.all<- maf.lt.all[,colnames(maf.lt.all)!="key"]
maf.lt.all[1:4,]

colnames(a.indel)[colnames(a.indel) == "MAF.lt:0.5"]<-"MAF.lt:0.05"
posns<-match(colnames(maf.lt.all),colnames(a.indel))
missing<-is.na(posns)
sum(missing)
## i<-1
## if(i in 1:sum(!missing)){
  
##  a.indel[, posns[!missing][i] ]  <- maf.lt.all[,!missing][i]
## }

dim(a.indel)
maf.lt.all[1:4,]
##AN: removing because it is transferred into "maf.lt."  i.e. :  maf.lt<-maf.lt & as.logical(target.table[,a.filter.cols.maf[imaff]])
if(exists("target.table")){rm(target.table)}


   
################### end get frequency table

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################


   





#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
length(all.possible.samples)
pheno.types

   
#####################################################
#################################### fix missiing gene names
#################################### got some missing gene names still.
## all.genes[grep("GTPBP4",names(all.genes))]
no.gene<-is.na(a.indel[,"Gene.Names"]) | a.indel[,"Gene.Names"]=="NA"
all.genes<-sort(table(a.indel[,"Gene.Names"]),decreasing=TRUE)

ens.to.hgnc<-read.delim(ensembleID.to.HGNC.map.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
ens.to.hgnc[1:5,]
genes<-unique(a.indel[no.gene,"Gene.Embl"])

posns<-match(genes,ens.to.hgnc[,"Ensembl.Gene.ID"])
missing<-is.na(posns)

ens.to.hgnc<-ens.to.hgnc[posns[!missing],]
dups<-duplicated(ens.to.hgnc[,"Ensembl.Gene.ID"])
ens.to.hgnc<-ens.to.hgnc[!dups,]
ens.to.hgnc[1:5,]

posns<-match(a.indel[no.gene,"Gene.Embl"],ens.to.hgnc[,"Ensembl.Gene.ID"])
missing<-is.na(posns)
sum(missing)

##AN: this is where mapped gene names from 'ens.to.hgnc' is mapped to missing gene nmames in a.indel
a.indel[no.gene,"Gene.Names"]<-ens.to.hgnc[posns,"HGNC.symbol"]

a.indel[no.gene,"Gene.Names"][1:10]
ens.to.hgnc[posns,"HGNC.symbol"][1:10]

no.gene<-is.na(a.indel[,"Gene.Names"]) | a.indel[,"Gene.Names"]=="NA"
a.indel[no.gene,"Gene.Names"]<-a.indel[no.gene,"Gene.Embl"]
a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
#a.indel[a.dash,][1:5,1:50]
a.indel[a.dash,"Gene.Names"]<-a.indel[a.dash,"Feature.Embl"]

a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
#a.indel[a.dash,][1:5,1:50]

final.names<-a.indel[a.dash,"knownGene::gene"]
final.names<-gsub("\\(dist=\\d*\\)","",final.names)
final.names<-gsub("\\(dist=NONE\\)","",final.names)
final.names[grepl("^NO_ANNOVAR",final.names)]<-"-"

a.indel[a.dash,"Gene.Names"]<-final.names
a.dash<-a.indel[,"Gene.Names"]=="-"
sum(a.dash)
a.indel[a.dash,c("chr","start")]
final.names<-build.key(a.indel[a.dash,],c("chr","start"))
a.indel[a.dash,"Gene.Names"]<-final.names

unannotated.hits<-a.dash

all.genes<-sort(table(a.indel[,"Gene.Names"]),decreasing=TRUE)
all.genes[1:30]
  ##    MUC4       TTN     MUC16     MUC5B     MUC12   FAM182B    ZNF717     NBPF1   FAM182A    AHNAK2 
  ##    1357       540       366       347       326       289       284       283       279       269 
  ##   MST1L     HYDIN       FLG     SYNE1    NBPF14 POM121L8P       NEB      MUC6 POM121L1P    MST1P2 
  ##     267       255       249       238       233       213       212       211       204       198 
  ## PDE4DIP  HLA-DRB1     CTBP2     MUC17     OBSCN  FLJ45445     USH2A     LRP1B      SSPO    DNAH11 
  ##     196       192       186       184       184       178       175       172       172       167
print("finished checking gene.names")
   

#grep("NOTCH1",names(all.genes))

common.hit.genes<-{}
## all.genes["SCN2A"]
####################################end  fix missiing gene names
###############################################

#########################################################################################################################
################################ make snpinfo file for association test
print("Build snpinfo")
   
snpinfo.raw<-cbind(key,a.indel[,"Gene.Names"],a.indel[,"Gene.Names"])
snpinfo.raw[1:5,]
tail(snpinfo.raw)
colnames(snpinfo.raw)<-c("Name","gene","cluster")

            poly.gene.site<-grep("::",snpinfo.raw[,"gene"])

            length(poly.gene.site) #### fix gene names that llok like gene1::gene2 - make 2 entries in snpinfo one for each
            ipoly<-2
            all.extra<-{}
            if(length(poly.gene.site)>0){
               
                for( ipoly in 1 :length(poly.gene.site)){
                    gene<-unlist(strsplit(snpinfo.raw[poly.gene.site[ipoly],"gene"],split="::") )
                    if(is.null(dim(all.extra))){
                    all.extra<-cbind(snpinfo.raw[poly.gene.site[ipoly],"Name"],gene,gene)
                    colnames(all.extra)<-colnames(snpinfo.raw)
                     }else{
                    extra<-cbind(snpinfo.raw[poly.gene.site[ipoly],"Name"],gene,gene)
                    colnames(extra)<-colnames(snpinfo.raw)
                    all.extra<-rbind(all.extra,extra)
                }
                       
 
            } }# has a poly

# all.extra[1:5,]

snpinfo.raw<-rbind(snpinfo.raw,all.extra)


colnames(snpinfo.raw)<-c("Name","gene","cluster")
dim(snpinfo.raw)
snpinfo.raw[1:5,]
dim(snpinfo.raw)


clusters<-read.delim(cluster.definition.file,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
clusters

#clusters.wanted<-c("Clinical","FANC - ACID")
clusters.wanted<-colnames(clusters)
ic<-1


gene.aliases<-read.delim(gene.symbol.file.for.clusters,header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE)
gene.aliases[1:5,]
gene.aliases<-unwind(gene.aliases, "Aliases",delimit=", ")
gene.aliases[1:5,]
###############  all.genes

   

## gene.aliases[grepl("SLX1",gene.aliases[,2]),]
############FIX MISSING gene names in pathways tests
#Check clusters have genes that are missing from the sequencing list then make sure have the approved gene name

   

ic<-5
recode<-{}
for(ic in 1:length(clusters.wanted)){
#  print(ic)
   cluster.genes<-clusters[,clusters.wanted[ic]]
   cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
   cluster.genes<-unique(cluster.genes)
#   all.genes[1:5]
   missing.name<-!(cluster.genes %in% names(all.genes))
   if(sum( missing.name)>0){
     posns<-match(cluster.genes[missing.name],gene.aliases[, "Aliases"])
     missing<-is.na(posns)
     

     
     if(sum(missing)>0){
       print(paste("in cluster",clusters.wanted[ic], " missing"))
       print(cluster.genes[missing.name][!missing])
     }
     recode<-cbind(cluster.genes[missing.name][!missing],gene.aliases[posns[!missing], "Approved.Symbol"])

     print(paste("recoded"))
     print(recode)

     
     colnames(recode)<-c("old","new")
 ###### transfer to new gene lists
     posns<-match(clusters[,clusters.wanted[ic]],recode[,"old"])
     missing<-is.na(posns)

    clusters[!missing,clusters.wanted[ic]]<-recode[posns[!missing],"new"] ### redefine the clusters
  }
  }
 #########################################################    
     

   
   


snpinfo<-snpinfo.raw
ic<-1
for(ic in 1:length(clusters.wanted)){
   cluster.genes<-clusters[,clusters.wanted[ic]]
   cluster.genes<-cluster.genes[cluster.genes!="" | is.na(cluster.genes)]
   print(paste(clusters.wanted[ic],paste(cluster.genes,collapse=","),sep=": "))
   extra<-snpinfo.raw[snpinfo.raw[,"gene"] %in%  cluster.genes,]
   if(!is.null(dim(extra))){
   print(dim(extra))
   extra[,"cluster"]<-clusters.wanted[ic]
   snpinfo<-rbind(snpinfo, extra)
  } ## extra contains some data

 }
snpinfo.ori<-snpinfo

   
#snpinfo.ori[1:20,]

### HAVE
# snpinfo.raw (original from a.indel)
# snpinfo # with extra clusters
# snpinfo.raw a permanent copy of snpinfo

## "FANCM " "MHF1"   "MHF2"   "FAAP24"
clusters[,1]
chk<-apply(clusters,2,function(x){ length(x[x!=""])})
write.table(clusters,file="/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/clusters_as.using.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

   
################################ make snpinfo file for association test
############################
#save(list=c("a.indel"),file="indels.RData")

   

############################################################################ do coverage filtered genotypes
   ############################################################################ do coverage filtered genotypes
num.cores
library(doMC)
registerDoMC(cores=num.cores)
## library("doParallel")
## registerDoParallel(cores=num.cores)
num.bits<-num.cores
#



the.samples<-colnames(a.indel)[grepl(".GT$", colnames(a.indel))]

while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out

# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8
## rownames(fil.genotypes)<-key
dim(fil.genotypes)
colnames(fil.genotypes)[1:5]
rownames(fil.genotypes)[1:5]

dim(fil.genotypes)
dim(a.indel)

tail(rownames(a.indel))

## a.indel<-a.indel[1:342791,]
## fil.genotypes<-fil.genotypes[1:342791,]

############################################################################ END do coverage filtered genotypes

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
############################### do one phenotype could be a loop #################### 
############################### do one phenotype ##################################################################
 ipheno<-1
# for(ipheno in 1:length(pheno.types)){  ## short circut look for now


###################################################################################################################################
##################################################### set association formula
print(paste("Doing phenotype:",pheno.types[ipheno]))
target.pheno<-names(pheno.types)[ipheno]
target.pheno.col<-pheno.types[ipheno]

formula<-paste(target.pheno.col,"~",paste(covars,collapse="+"),sep="")
print(formula)
formula<-formula(formula)

###################################################################################################################################
##################################################### END association formula
###############################################################




###################################################################################################################################
##################################################### ADD any extra cohort to pheno.ori that may want summary data for - summary.geno.extra

# sample.sheet.full<-sample.sheet.full[,1:31]
sample.sheet.full[1:2,]
colnames(sample.sheet.full)
sample.sheet.full[1:5,]
table(sample.sheet.full[,"SampleProject" ])

if ( !("SAMPLE" %in% colnames(sample.sheet.full))){
  SAMPLE<- sample.sheet.full[,"ParticipantCode"]
  sample.sheet.full<-cbind(sample.sheet.full,SAMPLE)
}

#if(sum(covars==1) & length(covars)==1){
 got.all.covars<-rep(TRUE,times=dim(sample.sheet.full)[1])
    
#}else{
#got.all.covars<-apply(sample.sheet.full[,covars],1,function(x) (sum(is.na(as.numeric(x))) ==0 ))
#}
all.possible.samples[(!(all.possible.samples %in% sample.sheet.full[,"SAMPLE"] ) )]
    

subset.samples<- (sample.sheet.full[,"SAMPLE"] %in% all.possible.samples ) & !is.na(sample.sheet.full[,target.pheno.col]) & sample.sheet.full[,target.pheno.col]!="NA"  & got.all.covars


    
sum(subset.samples)

    
##AN: pheno.or is the subset samples from sample sheet which contains different categories of samples

pheno.ori<-sample.sheet.full[subset.samples ,] ## pheno.ori only contains SAMPLES that have a phenotype
colnames(pheno.ori)
if(!("SAMPLE" %in% colnames(pheno.ori))){ ## add sample column if using PATIENT
  SAMPLE<-pheno.ori[,"PATIENT"]
  pheno.ori<-cbind(SAMPLE,pheno.ori)
 #  colnames(pheno.ori)[colnames(pheno.ori)=="PATIENT"]<- "SAMPLE"
}
colnames(pheno.ori)
print(dim(pheno.ori))
print(paste("Number Samples:",dim(pheno.ori)[1]))

dim(pheno.ori)
dim(sample.sheet.full)
length(all.possible.samples)

########## make sure sample sheet has no extras

table(pheno.ori[,"Project"])    


pheno.ori[1:5,]

###################################################################################################################################
##################################################### END  but phenotype file- pheno.ori from samplesheet for assoction
    




###################################################################################################################################
###################################################################################################################################
##################################################### END  but phenotype file- pheno.ori from samplesheet for assoction

    
the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")  ## samples same order as in pheno.ori
print(paste("Number samples: ",length(the.samples),sep=""))

#seq.type[1:57,]
posns<-match(pheno.ori[,"SAMPLE"],seq.type[,"Sample"])
missing<-is.na(posns)
sum(missing)
pheno.ori[missing,1:5]

if(!exists("pheno.ori")){
pheno.ori<-pheno
}
#pheno<-pheno.ori

#grep("AMLM12005R-G",seq.type[,"Sample"])

### Troels will update but all are Nextera
## capture<-seq.type[posns,"Capture.Method"]
## table(capture) ### all illume here
## pheno<-cbind(pheno,capture)
##       SampleProject FamilyCode SAMPLE PaternalID MaternalID
## NA.8        Control        ALL   3170          0          0
## NA.10       Control        ALL   3212          0          0
## NA.15       Control        ALL   3291          0          0
## NA.17       Control        ALL   3294          0          0
## NA.18       Control        ALL   3295          0          0
## NA.19       Control        ALL   3309          0          0

################################################################################

table(pheno.ori$SampleProject)

                         ## AML                    AML-Child       AML-NotDiagnosis-Child                    Asian-AML              Asian-AML-Child Asian-AML-NotDiagnosis-Child                Asian-Control                      Control 
                         ## 131                           17                           13                           15                            7                            5                           25                          429 

table(pheno.ori$Project)

    ## AMAS AOGC-NGS     MODY RSGB_AML     SDDS     SKDP TGCM-AML 
    ##   42      323       20       50       86       25       96
    
table(pheno.ori$capture)

## sum(!(pheno.ori.use[,"SAMPLE"]==pheno.ori[,"SAMPLE"]))
## the.projects<-c("cancer","Control","normal","PD","SCC","AK")

## cancer<-rep(FALSE,times=dim(pheno.ori)[1])
## cancer[pheno$SampleProject %in% c("AK","SCC")]<-TRUE

## normal<-rep(FALSE,times=dim(pheno.ori)[1])
## normal[pheno.ori$SampleProject %in% c("Normal")]<-TRUE

## Control<-rep(FALSE,times=dim(pheno.ori)[1])
## Control[pheno.ori$SampleProject %in% c("Control")]<-TRUE

## PD<-rep(FALSE,times=dim(pheno.ori)[1])
## PD[pheno.ori$SampleProject %in% c("PD")]<-TRUE

## SCC<-rep(FALSE,times=dim(pheno.ori)[1])
## SCC[grepl("_SCC$",pheno.ori$SAMPLE)]<-TRUE

## AK<-rep(FALSE,times=dim(pheno.ori)[1])
## AK[grepl("_AK1$",pheno.ori$SAMPLE) | grepl("_AK2$",pheno.ori$SAMPLE)]<-TRUE


## the.SCC<-pheno.ori$SAMPLE[grepl("_SCC$",pheno.ori$SAMPLE)]
## the.AK<-pheno.ori$SAMPLE[grepl("_AK1$",pheno.ori$SAMPLE) | grepl("_AK2$",pheno.ori$SAMPLE)]

## sum(PD)
## sum(normal)
## ####"LPH-001-27 Blood ok
## ## "LPH-001-27_PD"
## ## normal[pheno.ori$SAMPLE %in% c("LPH-001-27_PD")]<-TRUE
## ## PD[pheno.ori$SAMPLE %in% c("LPH-001-27_PD")]<-FALSE


## pheno.ori<-cbind(pheno.ori,cancer,Control,normal,PD,SCC,AK)

## pheno.ori[pheno.ori[,the.projects[1]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[2]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[3]],"SAMPLE"]
## pheno.ori[pheno.ori[,the.projects[4]],"SAMPLE"]

    
##NOVO
## [1] "cancer Num. samples: 74"
## [1] "Control Num. samples: 133"
## [1] "normal Num. samples: 25"
## [1] "PD Num. samples: 25"
## [1] "SCC Num. samples: 24"
## [1] "AK Num. samples: 50"

     ## AK Control  Normal      PD     SCC 
     ## 47     133      24      25      23 
#pheno.ori<-cbind(pheno.ori,SCC,AK)
############################################################## 
############## SET up groups to be analysed

########### special case if just using AOGC for QC
AOGC<-pheno.ori$Project=="AOGC-NGS" & !is.na(pheno.ori$Project)
pheno.ori$SAMPLE[AOGC]
#########################################################


the.projects<-sample.types
table(pheno.ori$SampleProject)
ist<-1

pheno.names.ori<-colnames(pheno.ori)
for(ist in 1:length(sample.types)){

a.type<-rep(FALSE,times=dim(pheno.ori)[1])
a.type[pheno.ori$SampleProject %in% sample.types[ist]]<-TRUE
assign(sample.types[ist],value=a.type)
pheno.ori<-cbind(pheno.ori,a.type)

     }
colnames(pheno.ori)<-c(pheno.names.ori,sample.types)
colnames(pheno.ori)


########### special case if just using AOGC for QC
sample.types.full<-sample.types
if(sum(AOGC)>0){
AOGC<-pheno.ori$Project=="AOGC-NGS" & !is.na(pheno.ori$Project)
pheno.ori$SAMPLE[AOGC]
pheno.ori<-cbind(pheno.ori,AOGC)
    sum(AOGC)
sample.types.full<-c(sample.types.full,"AOGC")
}


if("capture" %in% colnames(pheno.ori)){

the.captures<-table(pheno.ori$capture)
the.captures<-the.captures[the.captures>10]
icap<-1
for( icap in 1:length(the.captures)){
capture.case<-pheno.ori$capture==names(the.captures)[icap]  & pheno.ori[,cancer.group]
sum(capture.case)
pheno.ori<-cbind(pheno.ori,capture.case)
colnames(pheno.ori)[length(colnames(pheno.ori))]<-paste(names(the.captures)[icap],cancer.group,sep=".")
sample.types.full<-c(sample.types.full,paste(names(the.captures)[icap],cancer.group,sep="."))

capture.control<-pheno.ori$capture==names(the.captures)[icap]  & pheno.ori[,control.group]
sum(capture.control)
pheno.ori<-cbind(pheno.ori,capture.control)
colnames(pheno.ori)[length(colnames(pheno.ori))]<-paste(names(the.captures)[icap],control.group,sep=".")
sample.types.full<-c(sample.types.full,paste(names(the.captures)[icap],control.group,sep="."))

}

}    

#colnames(pheno.ori)[16:21]<-sample.types.full[4:9]

##AN: This is where Paul asked to add this PD
# PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]

## #################
## table(pheno.ori$capture)
## nextera.case<-pheno.ori$capture=="KapD:NimX3" | pheno.ori$capture=="TruD:NimX3" & pheno.ori$AML
## sum(nextera.case)
## pheno.ori<-cbind(pheno.ori,nextera.case)
## sample.types.full<-c(sample.types.full,"nextera.case")
## #################

## table(pheno.ori$capture)
## nextera.control<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$Control
## sum(nextera.control)
## pheno.ori<-cbind(pheno.ori,nextera.control)
## sample.types.full<-c(sample.types.full,"nextera.control")

## #################
## table(pheno.ori$capture)
## nextera.case<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$AML
## sum(nextera.case)
## pheno.ori<-cbind(pheno.ori,nextera.case)
## sample.types.full<-c(sample.types.full,"nextera.case")
## #################

## table(pheno.ori$capture)
## nextera.control<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$Control
## sum(nextera.control)
## pheno.ori<-cbind(pheno.ori,nextera.control)
## sample.types.full<-c(sample.types.full,"nextera.control")


## #################
## table(pheno.ori$capture)
## trueSeq.case<-pheno.ori$capture=="TruD:TruX" & pheno.ori$AML
## sum(trueSeq.case)
## pheno.ori<-cbind(pheno.ori,trueSeq.case)
## sample.types.full<-c(sample.types.full,"trueSeq.case")
## #################

## table(pheno.ori$capture)
## trueSeq.control<-pheno.ori$capture=="TruD:TruX" & pheno.ori$Control
## sum(trueSeq.control)
## pheno.ori<-cbind(pheno.ori,trueSeq.control)
## sample.types.full<-c(sample.types.full,"trueSeq.control")







####################3

###############################################################################
if("Aligner" %in% colnames(pheno.ori)){
table(pheno.ori$Aligner)
if(length(table(pheno.ori$Aligner))>1){
bwa.case<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$AML
sum(bwa.case)
pheno.ori<-cbind(pheno.ori,bwa.case)
sample.types.full<-c(sample.types.full,"bwa.case")

table(pheno.ori$Aligner)
novoalign.case<-pheno.ori$Aligner=="novoalign"  & pheno.ori$AML
sum(novoalign.case)
pheno.ori<-cbind(pheno.ori,novoalign.case)
sample.types.full<-c(sample.types.full,"novoalign.case")
}else{
  
  pheno.ori<-pheno.ori[,!(colnames(pheno.ori) %in% "Aligner")]  ##AN: the idea behind this aligner loop is if pheno.ori has more than two aligner, include bwa list, otherwise remove aligner column from pheno.ori. This line removes pheno.ori 
}
}

## table(pheno.ori$Aligner)
## bwa.control<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$Control
## sum(bwa.control) #0
## pheno.ori<-cbind(pheno.ori,bwa.case)
## sample.types.full<-c(sample.types,"bwa.case")
###############################################################################

###############################################################################

###############################################################################

## table(pheno.ori$AffectionStatus)
## sum(is.na((pheno.ori$AffectionStatus)))


#if("PD" %in% colnames(pheno.ori)){
    
PD<-pheno.ori$AffectionStatus==2 & !pheno.ori$AML
pheno.ori$SAMPLE[PD]
pheno.ori$SampleProject[PD]
pheno.ori<-cbind(pheno.ori,PD)
    sum(PD)
sample.types.full<-c(sample.types.full,"PD")



#################
##AN: labelling nextera and trueseq PD samples; ##Apr5 : Try changing this to NxtD:NxtXR
table(pheno.ori$capture)
nextera.PD<-pheno.ori$capture=="NxtD:NxtXR" & pheno.ori$PD
sum(nextera.PD)
pheno.ori<-cbind(pheno.ori,nextera.PD)
colnames(pheno.ori)[colnames(pheno.ori)=="nextera.PD"]<-"NxtD:NxtXR.PD"
sample.types.full<-c(sample.types.full,"NxtD:NxtXR.PD")
#################

#################
table(pheno.ori$capture)
trueSeq.PD<-pheno.ori$capture=="TruD:TruX" & pheno.ori$PD
sum(trueSeq.PD)
pheno.ori<-cbind(pheno.ori,trueSeq.PD)
colnames(pheno.ori)[colnames(pheno.ori)=="trueSeq.PD"]<-"TruD:TruX.PD"
sample.types.full<-c(sample.types.full,"TruD:TruX.PD")
#}
#################
##AN: labelling bwa and novoaligned PD samples

if("Aligner" %in% colnames(pheno.ori)){

table(pheno.ori$Aligner)
bwa.PD<-(pheno.ori$Aligner=="bwa" | pheno.ori$Aligner=="novoalign;bwa" ) & pheno.ori$PD
sum(bwa.PD)
pheno.ori<-cbind(pheno.ori,bwa.PD)
sample.types.full<-c(sample.types.full,"bwa.PD")


###############################################################################
table(pheno.ori$Aligner)
novoalign.PD<-pheno.ori$Aligner=="novoalign"  & pheno.ori$PD
sum(novoalign.PD)
pheno.ori<-cbind(pheno.ori,novoalign.PD)
sample.types.full<-c(sample.types.full,"novoalign.PD")
}



####################3

#########################################################



##############################################################
#439 unrelated germline controls ##############################################################
##############################################################
sample.types.full
the.projects.ori<-the.projects


the.projects<-sample.types.full
names(the.projects)<-the.projects
the.projects

    colnames(pheno.ori)
   
for (ir in 1: length(the.projects)){
  print(paste(the.projects[ir],"Num. samples:",sum(pheno.ori[,the.projects[ir]])))
      }

counts<-apply(pheno.ori[,the.projects],2,sum)
zero.samples<-names(counts)[counts<1]
the.projects<-the.projects[!(the.projects %in% zero.samples)]

# >counts
# [1] "Control Num. samples: 364"
# [1] "AML Num. samples: 131"
# [1] "AML-Child Num. samples: 17"
# [1] "AML-NotDiagnosis-Child Num. samples: 12"
# [1] "Asian-AML Num. samples: 15"
# [1] "Asian-AML-Child Num. samples: 7"
# [1] "Asian-AML-NotDiagnosis-Child Num. samples: 5"
# [1] "Asian-Control Num. samples: 6"
# [1] "AOGC Num. samples: 323"
# [1] "NxtD:NxtXR.AML Num. samples: 42"
# [1] "NxtD:NxtXR.Control Num. samples: 41"
# [1] "TruD:TruX.AML Num. samples: 89"
# [1] "TruD:TruX.Control Num. samples: 321"

## [1] "Control Num. samples: 429"
## [1] "AML Num. samples: 131"
## [1] "AML-Child Num. samples: 17"
## [1] "AML-NotDiagnosis-Child Num. samples: 13"
## [1] "Asian-AML Num. samples: 15"
## [1] "Asian-AML-Child Num. samples: 7"
## [1] "Asian-AML-NotDiagnosis-Child Num. samples: 5"
## [1] "Asian-Control Num. samples: 25"
## [1] "AOGC Num. samples: 323"


 pheno.ori[1:5,]

###################################################################################################################################
##################################################### END ADD any extra cohort to pheno.ori that may want summary data for - summary.geno.extra




###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
##################################################### BUILD summary.geno.extra
print("BUILD summary.geno.extra")
# summary.geno.extra.ori <-summary.geno.extra

summary.geno.extra<-{}
summary.geno.extra.ori<-{}
####################################################################################
#################################################################################### REGULAR
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
#targets<-sample.types # so include AOGC
targets



###### the.samples and pheno.ori in same order but the.samples has .GT extension.
it<-5

the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
print(targets[it])
#print(use.samples)
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it]," using ",length(use.samples)," samples ",sep=""))
} ## loop over targets
#################################################################################### FILTERED

#targets<-the.projects #c("NMD","ex.Control","AOGC")
targets
names(targets)<-paste(targets,".filt",sep="")
targets
it<-1
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
print(targets[it])
#print(use.samples)
length(use.samples)
genotypes<-fil.genotypes[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it]," using ",length(use.samples)," samples ",sep=""))
} ## loop over targets

#summary.geno.extra.ori<-summary.geno.extra
summary.geno.extra[1:5,]
colnames(summary.geno.extra)[grepl("GENO",colnames(summary.geno.extra))]


###################################################################################################################################
##################################################### END BUILD summary.geno.extra


###################################################################################################################################
#####################################################ADD genotype to a.indel.stats for later sunroutines

genotypes<-a.indel[,paste(all.possible.samples,"GT",sep=".")] #### use a.indel.ori otherwise
a.indel.stats<-cbind(a.indel.stats,genotypes)

## tail(colnames(a.indel.stats))
## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori


## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

######################################################################################################
######################################################################################################
##########################################################################
##########################################################################
########################### START GENOTYPE REVOVERY  MODEL ############################
##########################################################################
######################################################################################################
######################################################################################################


if(use.genotype.recovery & (!exists("a.indel.ori") | force.recovery.model)){ ### don't do revovery is alrready done!!
## cellularity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.cellularity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)

    
a.indel.ori<-a.indel
summary.geno.extra.ori<-summary.geno.extra


## alt.counts.thresh.ori<-alt.counts.thresh
## rare.in.Control.ori<-rare.in.Control
## rare.in.Control.filt.ori<-rare.in.Control.filt
## has.one.geno.ori<-has.one.geno
## rare.in.group.ori<-rare.in.group
## no.genotypes.filt.ori<-no.genotypes.filt
## no.genotypes.ori<-no.genotypes
## ok.missing.ori<-ok.missing


######################################################################################################
######################################################################################################
##################### Only revover rare genotypes likely to be somatic
n<-max(as.integer(summary.geno.extra.ori[,"TOTAL.Alleles.Control"]))

#p<-0.001
p<-0.01 ########### set MAF threshols HEREX1
sd.thresh<-6
n
p

alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh


## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

n<-as.numeric(summary.geno.extra.ori[,"TOTAL.Alleles.Control"])
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))
alt.counts.thresh.ori<-alt.counts.thresh


## alt.counts.thresh[1:50]
## summary.geno.extra.ori[1:5,grepl("GENO",colnames(summary.geno.extra))]
## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])< alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])< alt.counts.thresh

rare.in.Control<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.Control"])<= alt.counts.thresh
rare.in.Control.filt <-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh



sum(rare.in.Control)
sum(rare.in.Control.filt )
names(rare.in.Control)<-key
names(rare.in.Control.filt )<-key



## length(rare.in.Control)

## maf.lt.all[1:5,]
## maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])
maf.col<-paste("MAF.lt",p,sep=":")
maf.col
maf.filter<-as.logical(maf.lt.all[,maf.col])
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])
sum(maf.filter)
names(maf.filter)<-key
#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci



#pheno[pheno[,"AML"],c("SAMPLE","Blast.Count..", "Control","AML","SampleProject", "AffectionStatus"  )]
cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
# cellularity[pheno.ori[,"AML"],]
# cellularity[pheno.ori[,"Control"],]
# sample.types

PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"] | pheno.ori[,"Asian-Control"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")   table(pheno.ori[pheno.ori[,"AffectionStatus"]==2,"SAMPLE"] )
## PDs
#PDs.alt.counts<-alt.reads.reference.calls(a.indel,PDs,threshold=1)


cancer<-pheno.ori[pheno.ori[,"AML"],"SAMPLE"]
## cancer
#cancer.alt.counts<-alt.reads.reference.calls(a.indel,cancer,threshold=1)
#cancer.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,cancer,threshold=1)
Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]

#indels,the.samples,AD.extension="AD",threshold=1,prefix="",suffix=""

## a.indel.stats[1:5,1:20]
#Control.alt.counts<-alt.reads.reference.calls(a.indel.ori,Controls,AD.extension="AD",threshold=1)


Control.alt.counts<-alt.reads.reference.calls(a.indel.stats,Controls,AD.extension="FAD",threshold=1,prefix="",suffix="") ## needs a sample.GT column



## ##################################testing
## #Control.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,Control,threshold=1)
## #Control.alt.counts.ori<-Control.alt.counts
## Control.alt.counts1[1:5,]

## a.indel.ori[grepl("chr2:209113113:209113113",key),1:10]
## a.indel.ori[grepl("chr2:209113113:209113113",key),paste(cancer,".AD",sep="")]
## a.indel.ori[grepl("chr2:209113113:209113113",key),paste(cancer,".AD",sep="")]

## ## test<-c("chr12:11286221:11286221:T:A:snp","chr6:132029865:132029865:G:A:snp","chr7:82784807:82784807:C:A:snp","chr7:137206693:137206693:G:A:snp")

##  test<-c("chr15:90631934:90631934:C:T:snp","chr2:209113113:209113113:G:A:snp:209113113","chr15:90631838:90631838:C:T:snp") #,"chr12:11286221:11286221:T:A:snp")

##  test<-c("chr2:209113113:209113113:G:A:snp:209113113")
## test<-c("chr15:90631934:90631934:C:T:snp","chr15:90631838:90631838:C:T:snp")
## ## a.indel[test,paste(normals,".AD",sep="")]
## ## a.indel[test,paste(controls,".GT",sep="")]
## ## a.indel[test,paste(normals,".GT",sep="")]
## ## a.indel[test,paste(PDs,".AD",sep="")]
## ## a.indel[test,paste(PDs,".GT",sep="")]
## a.indel[test,paste(cancer,".AD",sep="")]
## a.indel[test,paste(cancer,".GT",sep="")]
## #a.indel.ori[test,paste(cancer,".GT",sep="")]

## #test<-c("chr15:40675107:40675107:C:T:snp") 
## #cbind(a.indel[test,paste(cancer,".GT",sep="")],a.indel.ori[test,paste(cancer,".GT",sep="")])
## Control.alt.counts[1:5,]
## geno.p1<-genotype.p.values(a.indel.ori[test,],c(PDs),AD.extension="AD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)

## geno.p0<-genotype.p.values(a.indel.ori[test,],c(cancer),AD.extension="AD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)
## geno.p1<-genotype.p.values(a.indel.stats[test,],c(cancer),AD.extension="FAD",Control.alt.counts[test,"Read.Balance"]/100,cellularity)
## #geno.p2<-genotype.p.values(a.indel[test,],c(PDs),normal.alt.counts[test,"Read.Balance"]/100,cellularity)


## 2*(pnorm(abs(c(4)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold=2*(pnorm(abs(c(4)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## p.threshold
## found.genotype<-  geno.p1<=p.threshold
## apply(found.genotype,1,sum)
## dim(geno.p1)
## apply(found.genotype,1,sum)


## rbind(a.indel.ori[test,paste0(c(PDs),".GT")][2,],
##       a.indel[test,paste0(c(PDs),".GT")][2,],
## a.indel[test,paste0(c(PDs),".AD")][2,],
##       signif(geno.p1[2,],2),
##       found.genotype[2,])


## geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
## geno.p[!found.genotype]<-"0/0"
## colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")
## colnames(genotypes)
## geno.p[test,]
## genotypes[test,colnames(geno.p)]



## use.samples<-paste(c(cancer),".GT",sep="")
## genotypes<-a.indel[test,use.samples]
## dim(genotypes)
## summary.geno<-genotype.summary(as.matrix(genotypes))
## colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"CANCER",sep=".")
## summary.geno



## ## use.samples<-paste(c(PDs),".GT",sep="")
## ## genotypes<-a.indel[test,use.samples]
## ## dim(genotypes)
## ## summary.geno<-genotype.summary(as.matrix(genotypes))
## ## colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),"PD",sep=".")
## ## summary.geno

## #summary.geno.extra[test,grepl("GENO",colnames(summary.geno.extra))]
## summary.geno.extra.ori[test,grepl("GENO",colnames(summary.geno.extra.ori))]


## ############## finish testing
## ##################################testing



# maf.filter
# rare.in.Control 
snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so does not work for SNPs
is.flat<-grepl("flat$",a.indel.ori[,"TYPE"])
has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0 | as.numeric(summary.geno.extra.ori[,"ALT.Alleles.PD"])>0
alt.count.thresh<-1  # use  p[true.sample.minor.counts <= alt.count.thresh]<-1 default more than one needed
alt.count.thresh.include<-2 # ignore above if  more that this number of samples have true calls (so if have 2 samples can recover samples with one call)
has.one.geno.ori<-has.one.geno
to.recover<-snp.only & has.one.geno & !is.flat # & maf.filter & rare.in.Control  ## maf.filter & rare.in.Control probably not required
recover.samples<-c(cancer,PDs,Controls)
sum(to.recover)
AD.lower.tail <-FALSE

# to.recover<- to.recover & a.indel[,"Gene.Names"] == "AMDHD2" 

#geno.p<-genotype.p.values(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail) ## COULD USE NORMAL HERE
geno.p<-genotype.p.values.row(a.indel.stats[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,alt.count.thresh.include,AD.lower.tail) ## COULD USE NORMAL HERE
                           ## Control.alt.counts[1:5,]
                           ##  normal.alt.counts[1:5,]


########################### TESTING AN example
## p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))


## dim(geno.p)
## length(pass)

## a.indel.sm<-a.indel[to.recover ,]
## ## test<-c("chr2:209113113:209113113:G:A:snp:209113113") # IDN
## ## test<-c("chrX:123195620:123195620:G:T:snp")#  $STAK2
## test<-a.indel.sm[,"Gene.Names"] == "AMDHD2" & pass[to.recover]

## a.indel.sm[test,1:5]
## p.threshold.z.thresh<-2
## positive<-geno.p[test,] <=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
## postive.samples<-apply(positive,2,function(x) sum(x,na.rm=TRUE) > 0)
## samples<-names(positive.samples)[positive.samples]
## #samples<-colnames(geno.p)[geno.p[test,] <=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))]

## samples
## geno.p[,samples]
## positive[,samples]
## sort(geno.p[4,])[1:20]

# samples<-c("80","45","AMLM12031CDF","AMLM12014N-R","88","AMLM12018AES","68")
## samples<-c("87","AMLM12003H-J","AMLM12003H-K","AMLM12015WPS","AMLM12021D-M","AMLM12027NM","AMLM12028DAK","AMLM12040A-J","AMLM12PAH037A-B")
## loc<-1
## trial<-rbind(a.indel[test,paste0(c(samples),".GT")],
##              a.indel.use[test,paste0(c(samples),".GT")],
##       a.indel[test,paste0(c(samples),".GQ")],
## a.indel.stats[test,paste0(c(samples),".FAD")],signif(geno.p[test,samples],2))
## trial
## chk<-grep("^GENO",colnames(summary.geno))
## trial

## range<-1:length(to.recover)
## recover.samples<-cancer
## a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] 
## a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][range,1:10] =="0/0"
## geno.p[80:100,1:10]
## a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][range,1:20] 
## a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][range,1:10] 

## range.c<-1:length(recover.samples) #30:50
## a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")][,range.c] 
## a.indel.use[to.recover ,paste(recover.samples,"GT",sep=".")][,range.c] =="0/0"
## geno.p[,range.c]
## a.indel.stats[to.recover ,paste(recover.samples,"FAD",sep=".")][,range.c] 
## a.indel.stats[to.recover ,paste(recover.samples,"DUP",sep=".")][,range.c]  


p.threshold.z.thresh<-4
p.threshold=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold #  0.000000001973175
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")


############# I don't want to correct GATK calls that are already 0/1 to 0/0:


genotypes<-a.indel.stats[to.recover ,paste(recover.samples,"GT",sep=".")] 
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]


sum(!( colnames(geno.p) %in% colnames(a.indel))) # must ge zero
to.transfer<-colnames(geno.p)[colnames(geno.p) %in% colnames(a.indel)]

posns<-match(rownames(a.indel),rownames(geno.p))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(geno.p)

a.indel[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]



## ########################### a.indel.stats contains non-recovered genotypes add the new ones now: would affest lib and aligner models slightly?
############################## problem is that bad aligner and lib position shoudl probbaly NOT be recovered so think this soudl remain as is
## ## genotypes<-a.indel[,paste(all.possible.samples,"GT",sep=".")] #### use a.indel.ori otherwise
## ## a.indel.stats<-cbind(a.indel.stats,genotypes)
## ## tail(colnames(a.indel.stats))

## posns<-match(rownames(a.indel.stats),rownames(geno.p))
## missing<-is.na(posns)
## sum(missing)
## sum(!missing)
## dim(geno.p)

## sum(!( colnames(geno.p) %in% colnames(a.indel.stats))) # must ge zero
## a.indel.stats[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]
## ########################################################




summary.geno.extra<-{}
####################################################################################
####################################################################################
#### MAY NEED TO ADJUST way use samples in selected based on selection below.
targets<-the.projects  #c("NMD","ex.Control","AOGC")
targets

###### the.samples and pheno in same order but the.samples has .GT extension.
it<-1

the.samples<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
print(targets[it])
print(use.samples)
length(use.samples)
genotypes<-a.indel[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it],sep=""))
} ## loop over targets
#################################################################################### FILTERED



################################################# redo filtered data


while((dim(a.indel)[1] %% num.bits)< 2){num.bits<-num.bits+1} ### go don't get matrix issues
num.bits
#(dim(a.indel)[1] %% num.bits)
## fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)

fil.genotypes<-foreach(a.indel.bit=iter(a.indel,by='row',chunksize=as.integer(dim(a.indel)[1]/num.bits) ), .combine='rbind', .multicombine=TRUE, .inorder=TRUE) %dopar% filtered.genotype(a.indel.bit,gsub(".GT$","",the.samples),prefix="",suffix="",20,0.02,0.98,0.20,0.80,7,2)
# 20,0.02,0.98,0.20,0.8,10,5 # for cancers where het may be way out

# 0.02< het(0/1)<0.98
# ref(0/0)< 0.2
# (1/1) > 0.8


rownames(fil.genotypes)<-key
dim(fil.genotypes)
## colnames(fil.genotypes)[1:5]
## rownames(fil.genotypes)[1:5]

dim(fil.genotypes)
dim(a.indel)

tail(rownames(a.indel))

###################


targets<-the.projects #c("NMD","ex.Control","AOGC")
targets
names(targets)<-paste(targets,".filt",sep="")
targets
it<-1
for(it in 1:length(targets)){
#use.samples<-the.samples[pheno.ori[,"SampleProject"]==targets[it]]
use.samples<-the.samples[pheno.ori[,targets[it]]]
print(targets[it])
print(use.samples)
length(use.samples)
genotypes<-fil.genotypes[,use.samples]
dim(genotypes)
summary.geno<-genotype.summary(as.matrix(genotypes))
colnames(summary.geno)<-paste(c("MAF","ALT.Alleles","REF.Alleles","TOTAL.Alleles","MISSING.Alleles","ALT_HOMO","ALT_HETRO","GENO"),names(targets)[it],sep=".")
#summary.geno[1:5,]
if(is.null(dim(summary.geno.extra))){
  summary.geno.extra<-summary.geno
}else{
  summary.geno.extra<-cbind(summary.geno.extra,summary.geno)
}
print(paste("Done: ",targets[it],sep=""))
} ## loop over targets



a.indel[1:5,1:10]
summary.geno.extra[1:5,]
colnames(summary.geno.extra)


rownames(summary.geno.extra)<-key

## ######################################### Check a few
## use.sample<-the.AK

## chk<-grep("^GENO",colnames(summary.geno))
## summary.geno["chr19:50169131:50169131:C:T:snp",chk]
## summary.geno.extra.ori<-summary.geno.extra
## a.snp<-"chr15:40675107:40675107:C:T:snp" # KNSTRN
## a.snp<-"chr19:50169131:50169131:C:T:snp" # BCL@L12
## a.snp<-"chr6:31940123:31940123:G:A:snp" #STK19

## chk<-grep("^GENO",colnames(summary.geno.extra))
## summary.geno.extra[a.snp,chk]

## ## chk<-grep("^GENO",colnames(summary.geno.extra.use))
## ## summary.geno.extra.use[a.snp,chk]

## chk<-grep("^GENO",colnames(summary.geno.extra.ori))
## summary.geno.extra.ori[a.snp,chk]







} #end genotype recovery

## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

############################################################################################
############################ END reCovery #############################
############################################################################################
##########################################################################################



############################################################################################
############################ START BWA NOVOALIGN FIX #############################
############################################################################################
##########################################################################################
#SNP_STATS:version1
#FILTER_SUMMARY:alt_mismatch5_proportion_gt_0.8;Strand_bias_proportion_gt_0.9;Read_end3_proportion_gt_0.8;High_not_called_alts_gte_3
#SUMMARY_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;prop_mismatch;prop_fstrand;prop_rstrand;prop_read_end
#SUMMARY_NOT_CALLED:ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count;high_alt_count
#SAMPLE.RF:GATK_called_snp;ref_count;alt_count;mismatch_alt_count;fstrand_alt_count;rstrand_alt_count;read_end_alt_count

## SDHA that align pooly with Nextera
## chr5:251680:251680:G:T:snp
## chr5:251672:251672:C:A:snp
## chr5:251674:251674:T:G:snp
## chr12:56676244:56676244:C:T:snp

## no support on ther strand for data
## chr5:251680:251680:G:T:snp
##
## test.samples<-c("AMLM12003H-J","AMLM12006CC","AMLM12007FF","AMLM12011J-G","AMLM12014DJG","AMLM12015WPS","AMLM12021D-M","AMLM12021M-R","AMLM12022N-A","AMLM12025GWR","AMLM12026MJC","AMLM12027NM","AMLM12028DAK","AMLM12030PGB","AMLM12031CDF","AMLM12036A-S","AMLM12036T-S","AMLM12038J-H","AMLM12040A-J","AMLM12PAH016K-B","AMLM12RMH026J-N")

## test.samples<-c("AMLM12001KP","AMLM12002K-B","AMLM12003H-J","AMLM12003H-K","AMLM12011J-G","AMLM12014N-R","AMLM12015WPS","AMLM12016D-F","AMLM12018AES","AMLM12018W-M","AMLM12019K-S","AMLM12019S-P","AMLM12020M-B","AMLM12021D-M","AMLM12021KS","AMLM12021M-R","AMLM12022N-A","AMLM12025GWR","AMLM12026MJC","AMLM12027MLM","AMLM12027NM","AMLM12028DAK","AMLM12029P-L","AMLM12030PGB","AMLM12031D-M","AMLM12032NRO","AMLM12036A-S","AMLM12036R-L","AMLM12036T-S","AMLM12037L-T","AMLM12037SAT","AMLM12038J-H","AMLM12039M-A","AMLM12040A-J","AMLM12PAH037A-B","AMLM12RMH026J-N")
## test.samples<-c("AMLM12001KP","AMLM12002K-B","AMLM12003H-J","100","35")
## test.samples<-c("AMLM12011J-G")

## test.loc<-c("chr5:251680:251680:G:T:snp","chr12:56676244:56676244:C:T:snp")
## get.samples<-expand.labels.to.samples(test.samples,c("FAD","TAD","DUP"),paste.after=TRUE)

## get.samples<-expand.labels.to.samples(test.samples,c("FAD"),paste.after=TRUE)
## get.samples<-expand.labels.to.samples(cancer,c("DUP"),paste.after=TRUE)

## rownames(filt)<-key
## filt[test.loc,1:12]

## grep(test.loc,key)

## a.indel.stats[test.loc,get.samples]            
## a.indel.stats[1:5,1:20]
######AN:Check if genotype recovery above needs to be used 

if("Aligner" %in% colnames(pheno.ori)){

bwa<-pheno.ori[pheno.ori[,"Aligner"] %in% c("bwa","novoalign;bwa")  ,"SAMPLE"] ## Leo Pharma
novoalign<-pheno.ori[ !(pheno.ori[,"Aligner"] %in% c("bwa","novoalign;bwa")) & pheno.ori[,"AffectionStatus"]!=1  ,"SAMPLE"]
#normals<-normals[normals!="LPH-001-27_PD"]

bwa
novoalign
#normal.alt.counts<-alt.reads.reference.calls(a.indel,novoalign,threshold=1)

novo.alt.counts<-alt.reads.reference.calls(a.indel.stats,novoalign,AD.extension="FAD",threshold=1)
bwa.alt.counts<-alt.reads.reference.calls(a.indel.stats,bwa,AD.extension="FAD",threshold=1) ## needs a sample.GT column


##  normal.allele.freq<-Control.alt.counts[test,"Read.Balance"]/100
## tumor.sample.major.counts<-normal.alt.counts[test,"REF.reads"]
## tumor.sample.minor.counts<-normal.alt.counts[test,"ALT.reads"]
## test<-c("chr19:2939268:2939289:ACCACCCTTACCCAAGGAGGCA:-:indel","chr13:21729956:21729956:A:G:snp","chr2:198363406:198363406:C:T:snp","chr22:20760282:20760282:A:C:snp")
#test<-c("chr5:251680:251680:G:T:snp","chr12:56676244:56676244:C:T:snp")


## poissonian.model.position(Control.alt.counts[test,"Read.Balance"]/100,normal.alt.counts[test,"REF.reads"], normal.alt.counts[test,"ALT.reads"],av.tumor.contamination=0)
## poss.model.ori<-poss.model
## pass.possonian.aligner.model.ori<-pass.possonian.aligner.model

poss.model<-poissonian.model.position(bwa.alt.counts[,"Read.Balance"]/100,novo.alt.counts[,"REF.reads"], novo.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## novo as reference negtive z more in bwa
## poss.model[test,]
poss.model[1:5,]
dim(poss.model)
dim(a.indel)
sum(rownames(bwa.alt.counts)!=rownames(novo.alt.counts))
sum(rownames(a.indel)!=rownames(poss.model))

#sum(poss.model[pass,"P-value"]<1e-5,na.rm=TRUE)
#poss.model[test,]  ### used in conjection with additional snp filtering
 # 2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) #2*(pnorm(abs(seq(from=3,to=5,by=0.1)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 2*(pnorm(abs(4.4), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)) ##4 sd greater than mean
#signif(2*(pnorm(-6, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)),digits=3)


pass.possonian.aligner.model<- as.numeric(poss.model[,"Z-score"]) >= -8  | is.na( poss.model[,"P-value"])

#poss.model[!pass.possonian.aligner.model,][1:10,]
#test %in% rownames(poss.model)[!pass.possonian.aligner.model]

length(pass.possonian.aligner.model)
#length(pass)
sum(!pass.possonian.aligner.model) #60 with subset 2
}else{ ## no aliger bias
pass.possonian.aligner.model<-rep(TRUE, dim(a.indel)[1])
}
    

################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test
################################# library test

#tapply(pheno.ori[,"capture"],pheno.ori[,"capture"],length)


colnames(pheno.ori)

## lib1<-pheno.ori[grepl("NxtXR",pheno.ori[,"capture"]) & pheno.ori[,"AffectionStatus"]!=1   ,"SAMPLE"] ## Leo Pharma
## are.lib1<-pheno.ori[,"SAMPLE"] %in% lib1

## lib2<-pheno.ori[ !(are.lib1) & pheno.ori[,"AffectionStatus"]!=1  ,"SAMPLE"]
## are.lib2<-pheno.ori[,"SAMPLE"] %in% lib2


################ TEST SAME CLASS DIFFERENT CAPTURE  ##AN: poissonian model to correct for capture bias for Control samples  -start
lib1<-pheno.ori[pheno.ori[,"NxtD:NxtXR.Control"],"SAMPLE"]
are.lib1<-pheno.ori[,"SAMPLE"] %in% lib1

lib2<-pheno.ori[pheno.ori[,"TruD:TruX.Control"],"SAMPLE"]
are.lib2<-pheno.ori[,"SAMPLE"] %in% lib2

table(pheno.ori[are.lib1,"capture"])
table(pheno.ori[are.lib1,"SampleProject"])
table(pheno.ori[are.lib2,"capture"])


##AN:check this funcrtion 'alt.reads.reference.calls' where it has: allele.depths.group[genotypes=="1/1" | genotypes=="0/1" | genotypes=="NA" | is.na(genotypes) ] <- NA
##AN:the idea is to compare homozygous reference genotypes only and if there is siginificant change between the capture groups. 
##AN:Read.Balance = (ALT.reads/(ALT.reads+REF.reads)) * 100
lib1.alt.counts<-alt.reads.reference.calls(a.indel.stats,lib1,AD.extension= FAD.col,threshold=1)
lib2.alt.counts<-alt.reads.reference.calls(a.indel.stats,lib2,AD.extension= FAD.col,threshold=1) ## needs a sample.GT column
#test<-c("chr7:44663912:44663912:A:T:snp","chrX:123195619:123195619:A:T:snp","chr7:44663913:44663913:C:T:snp","chrX:123200021:123200021:A:T:snp","chr5:251785:251785:G:A:snp")
## poss.model.lib.ori<-poss.model.lib
## pass.possonian.lib.model.ori<-pass.possonian.lib.model

    
poss.model.lib<-poissonian.model.position(lib2.alt.counts[,"Read.Balance"]/100,lib1.alt.counts[,"REF.reads"], lib1.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## lib1 had not ref reads if significant


#poss.model<-poissonian.model.position(lib1.alt.counts[,"Read.Balance"]/100,lib2.alt.counts[,"REF.reads"], lib2.alt.counts[,"ALT.reads"],1,lower.tail=FALSE) ## lib2 as reference negtive z more in bwa

#poss.model.lib[test,]


dim(poss.model.lib)
dim(a.indel.stats)

pass.possonian.lib.model<- (as.numeric(poss.model.lib[,"Z-score"]) <= 6  | is.na( poss.model.lib[,"P-value"]) ) | ( lib1.alt.counts[,"REF.reads"] < 100 | lib2.alt.counts[,"REF.reads"] < 100 )

#poss.model[!pass.possonian.aligner.model,][1:10,]
#test %in% rownames(poss.model)[!pass.possonian.aligner.model]
length(pass.possonian.lib.model)
#length(pass)
sum(!pass.possonian.lib.model) #739
#sum(!pass.possonian.lib.model.ori)


lib2.alt.counts[!pass.possonian.lib.model,][1:5,]

lib1.alt.counts[!pass.possonian.lib.model,][1:5,]


#summary.geno.extra[!pass.possonian.lib.model,grepl("GENO",colnames(summary.geno.extra))][1:5,]

##AN: poissonian model to correct for capture bias for Control samples  -end
##########################################################################
##########################################################################
############################ END BWA NOVOALIGN FIX #############################
##########################################################################
##########################################################################




################################################################################################################
################################################################################################################
################################ VARIANT dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################


## a.indel.use<-a.indel
## summary.geno.extra.use <-summary.geno.extra

## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori

## a.indel<-a.indel.use
## summary.geno.extra <-summary.geno.extra.use

###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols
###################################### GET HARDY Weinbery equilibrium of conrols

#colnames(summary.geno.extra)
#getHWE(obs_hets, obs_hom1, obs_hom2)
#hw.target<-"Control"  ## what to calculate HW with

hw.p.control<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,sep="")]) ## used 16 CPUs
hw.p.control.filt<-getHWE(summary.geno.extra[,paste("GENO.",hw.target,".filt",sep="")]) ## used 16 CPUs


length(hw.p.control)
names(hw.p.control)<-key
names(hw.p.control.filt)<-key

hw.p.control[1:5]
hw.p.control.filt[1:5]

hw.controls.ok<-hw.p.control > hwe.control.threshold
hw.controls.ok.filt<-hw.p.control.filt > hwe.control.threshold

######### testing
# 
# ######### testing
# hw.p.control[200100:200195]
# hw.p.control[22:25]
# summary.geno.extra[22:25,"GENO.Control"]
# test.loc<-c(24,200129:200143)
# test.loc<-c(1:10)
# geno<-strsplit(summary.geno.extra[test.loc,"GENO.Control"],split=",")
# AA<-as.numeric(unlist(lapply(geno,function(x) x[1])))
# BB<-as.numeric(unlist(lapply(geno,function(x) x[3])))
# AB<-as.numeric(unlist(lapply(geno,function(x) x[2])))
# x<-cbind(AA,AB,BB)
# x
# HWExactMat(x)
# hw.p.control[test.loc] ## correct in original order
# 
# ######
# 
# ######################## inbreeding and no genotypes ##AN: this is not used in the test
##AN: ask paul what these four lines below are
group.maf.thresh<-0.20 # group.maf.thresh<-0.10  if more common that this in Group then discard: Discard Rare but common in this cohort
missing.threshold<-0.20 # missing.threshold<-0.50  60 % of genotypes missing
missing.threshold.nimblgen<-0.50
missing.threshold.illumina<-0.50

targets<-the.projects

rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(targets.analysis,paste(targets.analysis,".filt",sep="")),sep="")]))
rare.in.group.table[1:5,]  ##AN: this what we get from groups in target.analysis 
# > rare.in.group.table[1:5,]  
# MAF.AML MAF.Control MAF.AML.filt MAF.Control.filt
# chr10:7830040:7830040:G:A:snp    TRUE        TRUE         TRUE             TRUE
# chr10:7830056:7830056:T:C:snp    TRUE        TRUE         TRUE             TRUE
# chr10:7830064:7830064:A:T:snp   FALSE       FALSE        FALSE             TRUE
# chr10:7830235:7830235:G:A:snp    TRUE        TRUE         TRUE             TRUE
# chr10:7830300:7830300:G:C:snp    TRUE        TRUE         TRUE             TRUE


# summary.geno.extra[1:5,]

## rare.in.group.table<-( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="") ]< group.maf.thresh) | ( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")] > (1-group.maf.thresh)) | (is.na( summary.geno.extra[,paste("MAF.",c(the.projects,paste(targets,".filt",sep="")),sep="")]))
## rare.in.group.table[1:5,]
## summary.geno.extra[1:5,]

#rare.in.group.test<-colnames(rare.in.group)
rare.in.group.test<-paste("MAF.",targets.analysis,sep="")
rare.in.group.test
##AN: ask Paul what rare.in.group is 
rare.in.group<-combine.boolean(rare.in.group.table,rare.in.group.test,"OR")
# > rare.in.group[1:4]
# chr10:7830040:7830040:G:A:snp chr10:7830056:7830056:T:C:snp chr10:7830064:7830064:A:T:snp chr10:7830235:7830235:G:A:snp 
# TRUE                          TRUE                         FALSE                          TRUE 

#sum(!rare.in.group)
#summary.geno.extra[!rare.in.group,paste("GENO.",targets.analysis,sep="")][1:5,]
#####AN: test for combine.boolean function
# function(a.matrix,columns.to.combine,combine.type,NA.set.to=TRUE){
#   ### Combine columns in a matrix using OR or AND per specified column OR USE cbind(vector1,vector2) to make a.matrix  with columns.to.combine of length one
#   ### Returns a bollean array
#   ############ set any NA's to true
#   if(length(columns.to.combine)==1){colnames(a.matrix)<-1:dim(a.matrix)[2];columns.to.combine<-colnames(a.matrix)}
#   print(columns.to.combine)
#   the.NAs<-is.na(a.matrix[,columns.to.combine])
#   count.col.NAs<-apply(as.matrix(the.NAs),2,sum)# one filter .col gives a vector
#   found.col.NAs<-count.col.NAs>0
#   if(sum(found.col.NAs)>0){
#     for(i in 1:length(found.col.NAs)){
#       if(found.col.NAs[i]){
#         a.matrix[,columns.to.combine[i]][is.na(a.matrix[,columns.to.combine[i]])]<-NA.set.to
#       }}}
#   ##########################
#   
#   
#   num.trues<-apply(as.matrix(a.matrix[,columns.to.combine]),1,function(x) sum(x,na.rm=TRUE))
#   if(combine.type=="OR"){the.filter<-num.trues >0} # ANY test is true the filter out 
#   if(combine.type=="AND"){the.filter<-num.trues==length(columns.to.combine)} # ALL test must be true the filter out 
#   the.filter
# }

rare.in.group[1:10]

####################
## rare.in.group<-combine.boolean(rare.in.group.table,c("MAF.AML","MAF.Control"),"OR")
## rare.in.group.filt<-combine.boolean(rare.in.group.table,c("MAF.AML.filt","MAF.Control.filt"),"OR")
## sum(!rare.in.group)
## sum(!rare.in.group.filt)

##AN: rare.in.group above is not used in the test END!!!!

##AN: no.genotypes is used for the analysis

no.genotypes.test<-paste("MAF.",targets.analysis,".filt",sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
summary.geno.extra[1:5,no.genotypes.test]
no.genotypes.filt<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes.filt[1:5]

no.genotypes.test<-paste("MAF.",targets.analysis,sep="")
no.genotypes.test
no.genotypes<-(summary.geno.extra[,no.genotypes.test]== 0)  | (is.na( summary.geno.extra[,no.genotypes.test])) # no genotypes in test classes for a mutataion after individaul quality filtering
no.genotypes[1:5,]
no.genotypes<-combine.boolean(no.genotypes,colnames(no.genotypes),"AND")
no.genotypes[1:5]



no.genotypes.extra.test<-paste("MAF.",cancer.group,sep="")
#no.genotypes.cancer<-(as.numeric(summary.geno.extra[,no.genotypes.cancer.test])== 0  & is.finite(as.numeric(summary.geno.extra[,no.genotypes.cancer.test]) ) )
no.genotypes.cancer<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.cancer)


if("PD" %in% colnames(pheno.ori)){

no.genotypes.extra.test<-paste("MAF.","PD",sep="")
no.genotypes.PD<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.PD)


no.genotypes.extra.test<-paste("MAF.","PD",".filt",sep="")
no.genotypes.PD.filt<-(summary.geno.extra[,no.genotypes.extra.test]== 0)               
sum(no.genotypes.PD.filt)



sum(no.genotypes)
    length(no.genotypes)

}else{ ### no PD class
    no.genotypes.PD<-rep(FALSE,times=dim(a.indel)[1])
    no.genotypes.PD.filt<-rep(FALSE,times=dim(a.indel)[1])
     }

#summary.geno.extra[1:5,]

print(paste0("High missing using: ",missing.targets))
high.missing<-{}
imt<-1
for(imt in 1:length(missing.targets)){
    the.missing.alleles<-paste("MISSING.Alleles",missing.targets[imt],sep=".")
    the.total.alleles<-paste("TOTAL.Alleles",missing.targets[imt],sep=".")
    ##AN: calculate missing alleles in cases and control
    a.missing.test<-as.numeric(summary.geno.extra[,the.missing.alleles])/(as.numeric(summary.geno.extra[,the.total.alleles])+as.numeric(summary.geno.extra[,the.missing.alleles]))
    
     if(is.null(length(high.missing)) | length(high.missing)==0){
       high.missing<-a.missing.test
       }else{
          high.missing<-cbind(high.missing,a.missing.test)  ##AN: here the first loop is "AML.filt" and second is "Control.filt"; see below it is renamed
        }
  }
high.missing.table<-high.missing
colnames(high.missing.table)<-missing.targets
rownames(high.missing.table)<-key

#   summary.geno.extra[1:5,c("GENO.AML.filt","GENO.Control","MISSING.Alleles.AML.filt","TOTAL.Alleles.AML.filt")]
#summary.geno.extra.ori[1:5,grep("AML$",colnames(summary.geno.extra.ori))]  
## high.missing<- cbind(as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH"])),
##                      as.numeric( summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.LOW.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.LOW.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.HIGH.pheno"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.HIGH.pheno"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.nimblegen"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.nimblegen"])),
##                      as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"])/(as.numeric(summary.geno.extra[,"TOTAL.Alleles.illumina"])+as.numeric(summary.geno.extra[,"MISSING.Alleles.illumina"]))
##                      )

colnames(high.missing.table)
high.missing.table[1:5,]

very.high.missing.controls.thresh<-0.75
very.high.missing.controls<- high.missing.table[,"Control.filt"]> very.high.missing.controls.thresh
sum(very.high.missing.controls)

high.total.missing<-subset(high.missing.table,select=missing.targets)
high.total.missing[1:5,]

high.total.missing<-high.total.missing <= missing.threshold
ok.missing<-combine.boolean(high.total.missing,missing.targets,"AND")
sum(ok.missing)

missing.targets
the.projects
colnames(high.missing.table)





################################################################################################################
################################################################################################################
################################ END VARIANT dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################





################################################################################################################
################################################################################################################
################################ START POSITION dependednt filters ######################################
################################################################################################################
################################################################################################################
################################################################################################################



#sum(high.total.missing | nimblegen.total.missing | illumina.total.missing)
#high.missing[1:5,]
################## get loci in common

###############################################################################################################################################
################# HOWEVER ALT ALLELES ALL IS USED AS MINOR ALLELE FREQUNCY  not alternative allele frequencies


# maf.lt.all<-a.indel[,colnames(a.indel)[grepl("MAF.lt:",colnames(a.indel))]]
#maf.lt.all[1:5,]
#maf.lt.all<-as.logical(maf.lt.all)
#as.logical(maf.lt.all[1:5,1])


####################################################################################
##################################### make the POSITION filter matrix QUALITY.FILTER FUNATIOAL group-MAF done here
   global.labs[!(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra))  )]
 if(sum( !(global.labs %in% c(colnames(a.indel),colnames(summary.geno.extra)))  )>0){print(paste("WARNING postion filters missing for data",
          toString(global.labs[!(global.labs %in% c(colnames(a.indel))) ])))}
#grep("ljb_gerp",colnames(a.indel))
#######AN:
# > global.labs
# [1] "QUAL"              "QD"                "HRun"              "SB"                "FILTER"            "PolyPhen.scores"   "SIFT.scores"       "mut.taster::score" "phylo::score"     
# [10] "PolyPhen.desc"     "SIFT.desc"         "GERP::score"       "MAF.ALL"           "MAF.HIGH"          "MAF.LOW"           "TYPE"             
#####
# > global.quality.cut
# QUAL                                     QD                                   HRun                                     SB 
# "50"                                  "0.5"                                    "5"                                    "1" 
# FILTER                                 FILTER                        PolyPhen.scores                        PolyPhen.scores 
# "PASS" "TruthSensitivityTranche99.90to100.00"                                  "0.1"                                  "0.4" 
# SIFT.scores                      mut.taster::score                           phylo::score                          PolyPhen.desc 
# "0.4"                                  "0.4"                                  "0.4"                             "damaging" 
# SIFT.desc                            GERP::score                            GERP::score                            GERP::score 
# "deleterious"                                  "2.5"                                    "2"                                    "0" 
# MAF.ALL                               MAF.HIGH                                MAF.LOW                                   TYPE 
# "0.25"                                 "0.25"                                 "0.25"                                 "flat" 
#####
# > global.quality.type
# QUAL              QD            HRun              SB     FILTER_PASS      FILTER_100    PolyPhen.low   PolyPhen.high       SIFT.high mut.taster.high      phylo.high    PolyPhen.bad 
# "numeric"       "numeric"       "numeric"       "numeric"        "factor"        "factor"       "numeric"       "numeric"       "numeric"       "numeric"       "numeric"        "factor" 
# SIFT.bad       GERP.high        GERP.low    GERP.unknown         MAF.ALL        MAF.HIGH         MAF.LOW            flat 
# "factor"       "numeric"       "numeric"       "numeric"       "numeric"       "numeric"       "numeric"        "factor" 

#[1] "WARNING postion filters missing for data HRun, phylo::score, MAF.HIGH, MAF.LOW"
############AN:END comment
##############$$$$$$$$$$$$$$ POSITION FILTER: combines a.indel, filter table and filter.table.pholy in this case.
qual<-position.quality.filter(  cbind( a.indel[,colnames(a.indel) %in% global.labs],summary.geno.extra[,colnames(summary.geno.extra) %in% global.labs] ), global.quality.cut,global.quality.type,global.quality.dirn)
dim(qual)
dim(a.indel)
rownames(qual)<-key
qual[1:5,]
#########$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


                      



the.types<-names(tapply(a.indel[,"ensGene::location"],a.indel[,"ensGene::location"],length))
the.types<-unique(unlist(strsplit(the.types,split=";")))
if( sum( !(the.types %in% c("NA","exonic",possible.mutations)))>0  ){print("WARNING ANNOVAR HAS NEW MUTATION TYPES DEFINED- REVIEW line 532")
                                                                print( the.types[!(the.types %in% c("exonic",possible.mutations))] )   }



#  filter.table.pholy[1:5,]
wanted.muts.coding<-test.for.coding.type(a.indel,geneanno.DB,interesting.coding.mutations)
sum(wanted.muts.coding)
   ####################### mhairi CHECK
wanted.muts.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],vep.coding,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## a.indel[wanted.muts.coding.vep,]
sum(wanted.muts.coding.vep)

  ## sum(wanted.muts.coding.vep)
  ##  sum(wanted.muts.coding)


###################### look for missense
missense.coding<-test.for.coding.type(a.indel,geneanno.DB,missense.variant)
sum(missense.coding)
   ####################### mhairi CHECK
missense.coding.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],missense.variant,delimit.by=",")  # filter.table.pholy[,"Consequence"] %in%  vep.coding
  ## a.indel[wanted.muts.coding.vep,]
sum(missense.coding.vep)

is.missense<-missense.coding.vep | missense.coding

## qual[1:50,"PolyPhen.low"] # true is above 0.1
## a.indel[1:50,"PolyPhen.scores"]
## a.indel[1:50,c("PolyPhen.desc","SIFT.desc")]
## qual[1:5,"SIFT.high"] 
## a.indel[1:5,"SIFT.scores"
synonymous<-test.for.coding.type(a.indel,geneanno.DB,synonymous.variant)
sum(synonymous)



is.benign.missense<-is.missense & !qual[,"PolyPhen.low"] & a.indel[,"PolyPhen.desc"] !="unknown" # chr9:21971016:21971016:G:A:snp

sum(!is.benign.missense) # is.benign.missense["chr9:21971016:21971016:G:A:snp"]
dim(a.indel)
###################### look for missense
###################### look for missense
###################### look for missense
##:AN
#PolyPhen-2 score
# The PolyPhen-2 score predicts the possible impact of an amino acid substitution on the structure and function of a human protein. This score represents the probability that a substitution is damaging. Ion Reporter™ Software reports the pph2-prob PolyPhen-2 score.
# The PolyPhen-2 score ranges from 0.0 (tolerated) to 1.0 (deleterious). Variants with scores of 0.0 are predicted to be benign. Values closer to 1.0 are more confidently predicted to be deleterious. The score can be interpreted as follows:
#   0.0 to 0.15 -- Variants with scores in this range are predicted to be benign.
# 0.15 to 1.0 -- Variants with scores in this range are possibly damaging.
# 0.85 to 1.0 -- Variants with scores in this range are more confidently predicted to be damaging.
# PolyPhen-2 and SIFT scores use the same range, 0.0 to 1.0, but with opposite meanings. A variant with a PolyPhen-2 score of 0.0 is predicted to be benign. A variant with a SIFT score of 1.0 is predicted to be benign.
##AN:end polyphen

############## noncoding double filtered for biotype ###
##AN: Interesting coding mutations
# > interesting.coding.mutations
# [1] "frameshift substitution"    "nonframeshift substitution" "nonframeshift deletion"     "nonframeshift insertion"    "frameshift deletion"        "frameshift insertion"      
# [7] "nonsynonymous SNV"          "stopgain SNV"               "stoploss SNV"               "splicing"                   "stopgain"                   "stoploss"  

interesting.NONcoding.mutations<-interesting.mutations[!(interesting.mutations %in% interesting.coding.mutations)] # "ncRNA_exonic"
wanted.muts.NONcoding<-test.for.coding.type(a.indel,geneanno.DB,interesting.NONcoding.mutations)
sum(wanted.muts.NONcoding)

wanted.interesting.to.prefilter<-test.for.coding.type(a.indel,geneanno.DB,interesting.to.prefilter) #"UTR3"      "UTR5"      "UTR5;UTR3" "snoRNA"
sum(wanted.interesting.to.prefilter)
# unique(a.indel[,"Consequence.Embl"])
  ####################### mhairi CHECK
wanted.interesting.to.prefilter.vep<-test.wanted.mutation(a.indel[,"Consequence.Embl"],vep.noncoding,delimit.by=",")  #filter.table.pholy[,"Consequence"] %in% vep.noncoding
sum(wanted.interesting.to.prefilter.vep)

wanted.muts.NONcoding.keep<-rep(FALSE,times=dim(a.indel)[1]) # "ncRNA_exonic" and "ncRNA_exonic" 
a.type<-wanted.noncoding.subtypes # list biotypes want to keep
for(itype in 1:length(a.type)){
  the.test<-grepl(a.type[itype],a.indel[,"Gene.Biotype"])
  the.test[is.na(the.test)]<-FALSE
  wanted.muts.NONcoding.keep<- wanted.muts.NONcoding.keep | the.test
}

 #  HERE wanted.muts.NONcoding.keep JUST DENOTES  "miRNA" &  "lincRNA" at this point USE wanted.muts.NONcoding to restrict to exones and splice BELOW
  wanted.muts.NONcoding.keep<-wanted.muts.NONcoding.keep & wanted.muts.NONcoding

dim(qual)
colnames(qual)
qual[1:5,]

bad.coding<-wanted.muts.coding | wanted.muts.coding.vep | (wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
sum(bad.coding)

bad.coding<-wanted.muts.coding | wanted.muts.coding.vep 
bad.non.coding<-(wanted.muts.NONcoding.keep & (qual[,"GERP.low"] | qual[,"GERP.unknown"]) ) | ( (wanted.interesting.to.prefilter | wanted.interesting.to.prefilter.vep) & (qual[,"GERP.high"]) )
bad.effect<-bad.coding | bad.non.coding

## table(a.indel[bad.non.coding,"Consequence.Embl"])
## test<-c("chr15:89848747:89848747:C:A:snp")
## test.pos<-grep(test,key)
## wanted.muts.NONcoding.keep[test.pos]
## wanted.interesting.to.prefilter[test.pos]
## wanted.interesting.to.prefilter.vep[test.pos]
## qual[test.pos,"GERP.low"]
## qual[test.pos,"GERP.unknown"]
## qual[test.pos,"GERP.high"]

# bad.coding<-test.for.coding.type(a.indel,geneanno.DB,c("stopgain SNV","stoploss SNV","frameshift deletion","frameshift insertion"))
#bad.frame<-test.for.coding.type(geneanno.table,geneanno.DB,c("frameshift deletion","frameshift insertion"))

## basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB"),"AND")
## gatk.qual<-combine.boolean(qual,c("FILTER_PASS", "FILTER_100" ),"OR")
## full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")

## basic.qual<-combine.boolean(qual,c("QUAL", "QD", "HRun", "SB","FILTER_100"),"AND")
## gatk.qual<-qual[,"FILTER_PASS"]
## full.qual<-combine.boolean(cbind(basic.qual,gatk.qual),"all","OR")

full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)

## full.qual<-gatk.qual
## sum(full.qual)
#full.qual<-qual[,"FILTER_PASS"]
sum(full.qual)



#any.functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","mut.taster.high","phylo.high","PolyPhen.bad","SIFT.bad","GERP.high","ljb_gerp.high"),"OR")
functional<-combine.boolean(qual,c("PolyPhen.low","SIFT.high","PolyPhen.bad","SIFT.bad","GERP.high"),"OR")
sum(functional)
functional<-functional | bad.coding # include really bad protein changes
sum(functional)

## chk<-"chr21:36252969:36252969:A:G:snp"
## a.indel[chk,1:50]

## maf.lt.all[1:5,]
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])

#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci

##########################################################################
##########################################################################



#################################################################


####################################
################Poly morphic SITE TESTS

REF.length<-nchar(as.character(a.indel[,"REF"]))
ALT.length<-nchar(as.character(a.indel[,"ALT"]))

large.indel<-REF.length>1 | ALT.length>1

are.repeats<-identify.repeats(a.indel,di.run.max=3,homo.run.max=5)

length(large.indel) 
length(are.repeats)
sum(are.repeats)
#rownames(a.indel)[are.repeats][1:20]

#################### in repeats looking  forward

#chk.in.repeat<-large.indel & !are.repeats
chk.in.repeat<- !are.repeats

are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="forward",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.forward<- key %in% remove.repeats

#remove.repeats[1:20]
sum(are.in.repeats.forward)
## [1] 6988

###################### in repeats looking back are.repeats[789]

sum(chk.in.repeat)
#chk.in.repeat<-large.indel & !are.repeats & !are.in.repeats.forward
chk.in.repeat<- !are.repeats & !are.in.repeats.forward

are.sub.repeat<-indentify.IN.repeat(a.indel[chk.in.repeat,],looking="back",bases.about=6,di.run.max=3,homo.run.max=5,genome="BSgenome.Hsapiens.UCSC.hg19")
remove.repeats<-key[chk.in.repeat][are.sub.repeat]
are.in.repeats.back<- key %in% remove.repeats

#remove.repeats[1:20]
sum(are.in.repeats.back)
## [1] 3224

are.in.repeats<- are.in.repeats.back | are.in.repeats.forward

length(are.in.repeats)
sum(are.in.repeats)



##########################################################################
##########################################################################
##########################################################################
not.flat.genotype<-!qual[,"flat"]
sum(not.flat.genotype)

colnames(qual)

is.unwound.geno<-grepl("snp:\\d+$",a.indel[,"TYPE"]) | grepl("indel:\\d+$",a.indel[,"TYPE"])



#hw.controls.ok[loci]

in.common.hit.gene <- a.indel[,"Gene.Names"] %in% common.hit.genes
in.common.hit.gene[1:5] 
sum(in.common.hit.gene)
the.chr
on.x.y<-a.indel[,"chr"] %in% c("X","Y","23","24","chrX","chrY")
sum(on.x.y)

#table(a.indel[,"TYPE"])
snp.only<-grepl("^snp",a.indel[,"TYPE"])



########################################################################## /media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/Gerp.test.r

## in.any.normal<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal"]) > 0
## in.any.normal.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal.filt"]) > 0

## summary.geno.extra[!in.any.normal.filt,grepl("normal",colnames(summary.geno.extra))][1:5,]
## sum(!in.any.normal)

## [1] "Control Num. samples: 173"
## [1] "cancer Num. samples: 14"
## [1] "normal Num. samples: 14"
## [1] "KapD:NimX3.cancer Num. samples: 14"
## [1] "KapD:NimX3.Control Num. samples: 52"
## [1] "TruD:NimX3.cancer Num. samples: 0"
## [1] "TruD:NimX3.Control Num. samples: 78"
## [1] "TruD:TruX.cancer Num. samples: 0"

icap<-1
bias.cancer<-{}
bias.labels<-{}
bias.cancer.pval<-{}
for( icap in 1:length(the.captures)){
  a.capture<-the.captures[icap]
  
  if( !(paste("ALT.Alleles",names(a.capture),cancer.group,sep=".") %in% colnames(summary.geno.extra)) ) {print("skip"); print(icap) ;next} ## no samples in that type
  
  
  ####################################### nextera bias
  n.capture.cancer<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles",names(a.capture),cancer.group,sep=".")]))
  n.cancer<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles.",cancer.group,sep="")]))
  
  p.default<-n.capture.cancer/n.cancer
  
  p.cohort<-as.integer(summary.geno.extra[,paste("TOTAL.Alleles",names(a.capture),cancer.group,sep=".")])/as.integer(summary.geno.extra[,paste("TOTAL.Alleles.",cancer.group,sep="")])
  p.cohort[!is.finite(p.cohort)]<-p.default
  names(p.cohort)<-rownames(summary.geno.extra)
  test<-dbinom(as.integer(summary.geno.extra[,paste("ALT.Alleles",names(a.capture),cancer.group,sep=".")]), as.integer(summary.geno.extra[,paste("ALT.Alleles.",cancer.group,sep="")]), p.cohort, log=FALSE)
  #test.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## cancer where ther are zero trials
  names(test)<-rownames(summary.geno.extra)
  a.bias<-(as.numeric(test) < 0.001)
  sum(a.bias)
  if( is.null((bias.cancer)) ){bias.cancer<-a.bias; bias.labels<-names(a.capture);  print(icap) ;bias.cancer.pval<-test  }else{
    bias.cancer<-cbind(bias.cancer,a.bias); bias.labels<-cbind(bias.labels,names(a.capture)) ; bias.cancer.pval<-cbind(bias.cancer.pval,test) }
  #######################################################
  
}

if(!is.null(dim(bias.cancer))){
  colnames(bias.cancer)<-paste(bias.labels,cancer.group,sep=".")
  colnames(bias.cancer.pval)<-paste(bias.labels,cancer.group,sep=".")
  rownames(bias.cancer)<-key
  rownames(bias.cancer.pval)<-key
  fail.bias.cancer<-combine.boolean(bias.cancer,colnames(bias.cancer),"AND")
}else{
  # bias.cancer[1:5]
  bias.cancer.pval<-as.matrix(bias.cancer.pval)
  colnames(bias.cancer.pval)<-paste(bias.labels,cancer.group,sep=".")
  rownames(bias.cancer.pval)<-key
  fail.bias.cancer<-bias.cancer
}


########################### Controls

icap<-1
bias.control<-{}
bias.labels<-{}
bias.control.pval<-{}
for( icap in 1:length(the.captures)){
  a.capture<-the.captures[icap]
  
  if( !(paste("ALT.Alleles",names(a.capture),control.group,sep=".") %in% colnames(summary.geno.extra))) {print("skip"); print(icap) ;next} ## no samples in that type
  
  
  ####################################### nextera bias
  n.capture.control<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles",names(a.capture),control.group,sep=".") ]))
  n.control<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles.",control.group,sep="")]))
  
  p.default<-n.capture.control/n.control
  
  p.cohort<-as.integer(summary.geno.extra[,paste("TOTAL.Alleles",names(a.capture),control.group,sep=".")]) / as.integer(summary.geno.extra[,paste("TOTAL.Alleles.",control.group,sep="")])
  p.cohort[!is.finite(p.cohort)]<-p.default
  names(p.cohort)<-rownames(summary.geno.extra)
  test<-dbinom(as.integer(summary.geno.extra[,paste("ALT.Alleles",names(a.capture),control.group,sep=".")]), as.integer(summary.geno.extra[,paste("ALT.Alleles.",control.group,sep="")]), p.cohort, log=FALSE)
  #test.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## control where ther are zero trials
  names(test)<-rownames(summary.geno.extra)
  a.bias<-(as.numeric(test) < 0.001)
  sum(a.bias)
  if( is.null((bias.control)) ){bias.control<-a.bias; bias.labels<-names(a.capture); print(icap) ; bias.control.pval<-test }else{
    bias.control<-cbind(bias.control,a.bias) ; bias.labels<-cbind(bias.labels,names(a.capture)) ; bias.control.pval<-cbind(bias.control.pval,test) }
  #######################################################
  
}

if(!is.null(dim(bias.control))){
  bias.control[1:5,]
  colnames(bias.control)<-paste(bias.labels,control.group,sep=".")
  colnames(bias.control.pval)<-paste(bias.labels,control.group,sep=".")
  rownames(bias.control)<-key
  rownames(bias.control.pval)<-key
  fail.bias.control<-combine.boolean(bias.control,colnames(bias.control),"AND")
}else{
  bias.control[1:5]
  bias.control.pval<-as.matrix(bias.control.pval)
  colnames(bias.control.pval)<-paste(bias.labels,control.group,sep=".")
  rownames(bias.control.pval)<-key
  fail.bias.control<-bias.control
  
}


sum(fail.bias.control)

############## PDs

# [1] "skip"
# [1] 1
# [1] "skip"
# [1] 2
# Error in as.matrix(bias.PD.pval) : 
#   error in evaluating the argument 'x' in selecting a method for function 'as.matrix': Error: object 'bias.PD.pval' not found


if("PD" %in% colnames(pheno.ori)){
  icap<-1
  bias.PD<-{}
  bias.labels<-{}
  bias.labels.pval<-{}
  for( icap in 1:length(the.captures)){
    a.capture<-the.captures[icap]
    
    if( !(paste("ALT.Alleles.",names(a.capture),".PD",sep="") %in% colnames(summary.geno.extra))) {print("skip"); print(icap) ;next} ## no samples in that type
    
    
    ####################################### nextera bias
    n.capture.PD<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles.",names(a.capture),".PD",sep="") ]))
    n.PD<-max(as.integer(summary.geno.extra[,paste("ALT.Alleles.",PD.group,sep="")]))
    
    p.default<-n.capture.PD/n.PD
    
    p.cohort<-as.integer(summary.geno.extra[,paste("TOTAL.Alleles.",names(a.capture),".PD",sep="")]) / as.integer(summary.geno.extra[,paste("TOTAL.Alleles.",PD.group,sep="")])
    p.cohort[!is.finite(p.cohort)]<-p.default
    names(p.cohort)<-rownames(summary.geno.extra)
    test<-dbinom(as.integer(summary.geno.extra[,paste("ALT.Alleles.",names(a.capture),".PD",sep="")]), as.integer(summary.geno.extra[,paste("ALT.Alleles.",PD.group,sep="")]), p.cohort, log=FALSE)
    #test.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## PD where ther are zero trials
    names(test)<-rownames(summary.geno.extra)
    a.bias<-(as.numeric(test) < 0.001)
    sum(a.bias)
    if( is.null((bias.PD)) ){bias.PD<-a.bias; bias.labels<-names(a.capture); print(icap) ; bias.PD.pval<-test  }else{
      bias.PD<-cbind(bias.PD,a.bias) ; bias.labels<-cbind(bias.labels,names(a.capture)) ; bias.PD.pval<-cbind(bias.PD.pval,test) }
    #######################################################
    
  }
  
  if(!is.null(dim(bias.PD))){
    bias.PD[1:5,]
    colnames(bias.PD)<-paste(bias.labels,PD.group,sep=".")
    colnames(bias.PD.pval)<-paste(bias.labels,PD.group,sep=".")
    rownames(bias.PD)<-key
    rownames(bias.PD.pval)<-key
    fail.bias.PD<-combine.boolean(bias.PD,colnames(bias.PD),"AND")
  }else{
    bias.PD[1:5]
    bias.PD.pval<-as.matrix(bias.PD.pval)
    colnames(bias.PD.pval)<-paste(bias.labels,PD.group,sep=".")
    rownames(bias.PD.pval)<-key
    fail.bias.PD<-bias.PD
    
  }
  
  sum(fail.bias.PD)
  
}else{ # no PD
  fail.bias.PD<-rep(FALSE,times=dim(a.indel)[1])
  bias.PD<-rep(FALSE,times=dim(a.indel)[1])
  bias.PD.pval<-rep(1,times=dim(a.indel)[1])
  bias.PD.pval<-as.matrix(bias.PD.pval)
  colnames(bias.PD.pval)<-paste("none",PD.group,sep=".")
  rownames(bias.PD.pval)<-key
  fail.bias.PD<-bias.PD
}

## ####################################### nextera bias
## n.nextera.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.case"]))
## n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

## p.nextera.default<-n.nextera.case/n.case

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.nextera.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
## p.cohort[!is.finite(p.cohort)]<-p.nextera.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.nextera<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
## #test.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
## names(test.nextera)<-rownames(summary.geno.extra)
## nextera.bias<-(as.numeric(test.nextera) < 0.001)
## #######################################################





## ####################################### trueSeq bias
## n.trueSeq.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.case"]))
## n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

## p.trueSeq.default<-n.trueSeq.case/n.case

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.trueSeq.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
## p.cohort[!is.finite(p.cohort)]<-p.trueSeq.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.trueSeq<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
## #test.trueSeq[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
## names(test.trueSeq)<-rownames(summary.geno.extra)
## trueSeq.bias<-(as.numeric(test.trueSeq) < 0.001)
## #######################################################


## ####################################### bwa bias
## n.bwa.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.case"]))
## n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

## p.bwa.default<-n.bwa.case/n.case

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.bwa.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
## p.cohort[!is.finite(p.cohort)]<-p.bwa.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.bwa<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
## #test.bwa[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
## names(test.bwa)<-rownames(summary.geno.extra)
## bwa.bias<-(as.numeric(test.bwa) < 0.001)
## #######################################################

## ####################################### novoalign bias
## n.novoalign.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.case"]))
## n.case<-max(as.integer(summary.geno.extra[,"ALT.Alleles.AML"]))

## p.novoalign.default<-n.novoalign.case/n.case

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.novoalign.case"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.AML"])
## p.cohort[!is.finite(p.cohort)]<-p.novoalign.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.novoalign<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.case"]), as.integer(summary.geno.extra[,"ALT.Alleles.AML"]), p.cohort, log=FALSE)
## #test.novoalign[as.integer(summary.geno.extra[,"ALT.Alleles.AML"])==0]<-1 ## case where ther are zero trials
## names(test.novoalign)<-rownames(summary.geno.extra)
## novoalign.bias<-(as.numeric(test.novoalign) < 0.001)
## #######################################################




## ############################################### repeat for the PD cohort
## ####################################### nextera bias
## n.nextera.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.PD"]))
## n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

## p.nextera.default<-n.nextera.PD/n.PD

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.nextera.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
## p.cohort[!is.finite(p.cohort)]<-p.nextera.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.PD.nextera<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.nextera.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
## #test.PD.nextera[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
## names(test.PD.nextera)<-rownames(summary.geno.extra)
## nextera.PD.bias<-(as.numeric(test.PD.nextera) < 0.001)
## #######################################################

## ####################################### trueSeq bias
## n.trueSeq.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.PD"]))
## n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

## p.trueSeq.default<-n.trueSeq.PD/n.PD

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.trueSeq.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
## p.cohort[!is.finite(p.cohort)]<-p.trueSeq.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.PD.trueSeq<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.trueSeq.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
## #test.PD.trueSeq[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
## names(test.PD.trueSeq)<-rownames(summary.geno.extra)
## trueSeq.PD.bias<-(as.numeric(test.PD.trueSeq) < 0.001)
## #######################################################


## ####################################### bwa bias
## n.bwa.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.PD"]))
## n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

## p.bwa.default<-n.bwa.PD/n.PD

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.bwa.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
## p.cohort[!is.finite(p.cohort)]<-p.bwa.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.PD.bwa<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.bwa.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
## #test.PD.bwa[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
## names(test.PD.bwa)<-rownames(summary.geno.extra)
## bwa.PD.bias<-(as.numeric(test.PD.bwa) < 0.001)
## #######################################################

## ####################################### novoalign bias
## n.novoalign.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.PD"]))
## n.PD<-max(as.integer(summary.geno.extra[,"ALT.Alleles.PD"]))

## p.novoalign.default<-n.novoalign.PD/n.PD

## p.cohort<-as.integer(summary.geno.extra[,"TOTAL.Alleles.novoalign.PD"])/as.integer(summary.geno.extra[,"TOTAL.Alleles.PD"])
## p.cohort[!is.finite(p.cohort)]<-p.novoalign.default
## names(p.cohort)<-rownames(summary.geno.extra)
## test.PD.novoalign<-dbinom(as.integer(summary.geno.extra[,"ALT.Alleles.novoalign.PD"]), as.integer(summary.geno.extra[,"ALT.Alleles.PD"]), p.cohort, log=FALSE)
## #test.PD.novoalign[as.integer(summary.geno.extra[,"ALT.Alleles.PD"])==0]<-1 ## case where ther are zero trials
## names(test.PD.novoalign)<-rownames(summary.geno.extra)
## novoalign.PD.bias<-(as.numeric(test.PD.novoalign) < 0.001)
## #######################################################



use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,cancer.group]]
the.QG.cancer<-genotype.GQ.summary(a.indel,use.samples.GQ,cancer.group)

use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,control.group]]
the.QG.Controls<-genotype.GQ.summary(a.indel,use.samples.GQ,control.group)



####
if("PD" %in% colnames(pheno.ori)){
use.samples.GQ<-paste(pheno.ori[,"SAMPLE"],"GT",sep=".")
use.samples.GQ<-use.samples.GQ[pheno.ori[,"PD"]]
the.QG.PD<-genotype.GQ.summary(a.indel,use.samples.GQ,"PD")
}else{  ### no PD
    the.QG.PD<-matrix(data=NA,nrow=dim(the.QG.cancer)[1],ncol=dim(the.QG.cancer)[2])
    colnames(the.QG.PD)<-gsub(".cancer$",".PD",colnames(the.QG.cancer))
}
    


GQ.cancer.pass<-the.QG.cancer[,paste("the.mean.",cancer.group,sep="")] > 50 | is.na(the.QG.cancer[,paste("the.mean.",cancer.group,sep="")])
GQ.Control.pass<-the.QG.Controls[,paste("the.mean.",control.group,sep="")] > 50 | is.na(the.QG.Controls[,paste("the.mean.",control.group,sep="")])

##AN: commented out this test for PD, should be the.GQ.PD
# if("PD" %in% colnames(pheno.ori)){
# GQ.PD.pass<-the.QG.Controls[,paste("the.mean.","PD",sep="")] > 50 | is.na(the.QG.Controls[,paste("the.mean.","PD",sep="")])
# }else{  ### no PD
#     GQ.PD.pass<-rep(TRUE,times=dim(a.indel)[1])
#      names(GQ.PD.pass)<-key
# }


if("PD" %in% colnames(pheno.ori)){
  GQ.PD.pass<-the.QG.PD[,paste("the.mean.","PD",sep="")] > 50 | is.na(the.QG.PD[,paste("the.mean.","PD",sep="")])
}else{  ### no PD
  GQ.PD.pass<-rep(TRUE,times=dim(a.indel)[1])
  names(GQ.PD.pass)<-key
}




##AN: this part is not used in the test END



#########################  Choose Somatic
###AN: Ask paul what this here is

if(subset.somatic){

in.any.normal<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal"]) > 0
in.any.normal.filt<-as.numeric(summary.geno.extra[,"ALT.Alleles.normal.filt"]) > 0
names(in.any.normal)<-key
names(in.any.normal.filt)<-key

## in.any.normal.ori<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.normal"]) > 0
## in.any.normal.filt<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.normal.filt"]) > 0
#summary.geno.extra[!in.any.normal.filt,grepl("normal",colnames(summary.geno.extra))][1:5,]
sum(!in.any.normal)
sum(in.any.normal.filt & !in.any.normal) # 0 must be zero

}else{
    in.any.normal<-rep(FALSE,times=dim(a.indel)[1])
    names(in.any.normal)<-key
     in.any.normal.filt<-rep(FALSE,times=dim(a.indel)[1])
     names(in.any.normal.filt)<-key
}
    
##########################################################################
##########################################################################
##########################################################################
icc<-3
 if(target.pheno.col %in% case.control){
   for(icc in 1:length(case.control.classes)){
     recode<-  pheno.ori[,target.pheno.col] %in% names(case.control.classes)[icc]
     pheno.ori[recode,target.pheno.col]<-as.numeric(case.control.classes[icc])
   }}
  
pheno.ori[,target.pheno.col]<-as.numeric(pheno.ori[,target.pheno.col])
formula



## pheno.use<-pheno[pheno[,target.pheno] %in% c(0,1),]
## dim(pheno.use) pheno.ori.ori<-pheno.ori
if(!exists("pheno.ori")){
pheno.ori<-pheno
}

## if(!identical(colnames(pheno),colnames(pheno.ori))){
## pheno.ori<-pheno
## }

################
################
  
pheno<-pheno.ori[pheno.ori[, "SampleProject"] %in% c(0,1) ,] #SCC dim(pheno.ori)
 table(pheno[, "SampleProject"])

pheno.ori[,c("SAMPLE","Project","SampleProject")]

table(pheno.ori[, "Project"])
pheno<-pheno.ori[ pheno.ori[, "SampleProject"] %in% c(0,1) & !(pheno.ori[, "Project"] %in% cohorts.to.exclude.from.assoc )  ,] # & !is.na((pheno.ori[, "Project"])) ### NA in project pass through
table(pheno[, "SampleProject"])
table(pheno[, "Project"])
#pheno<-pheno.ori
## pheno[!(pheno[, "Project"] %in% c("MODY","SDDS","SKDP"),"SAMPLE"]
##  pheno.ori(pheno.ori[, "Project"] %in% c("AOGC-NGS")) ,c("SAMPLE","SampleProject")]
## pheno[pheno[, "SampleProject"] %in% c(0),c("SAMPLE","SampleProject","Project")]

##AN: pheno.or[,"SampleProject"] is where the case control is classified
##AN:saving pheno for another test

the.samples.use<-pheno[,"SAMPLE"]
the.samples.use<-paste(the.samples.use,".GT",sep="")






## heno[pheno[, "capture"]=="NA",]

 table(pheno[, "SampleProject"])
##   0   1 
## 323 135 
## table(pheno[, "capture"])
## KapD:NimX3           NxtD:NxtX            NxtD:NxtXR    TruD:NimX     TruD:NimX3;NxtD:NxtXR     TruD:TruX  Unknown
##     1                  0                     161            3                     1                183         6

table(pheno[, "Project"])
## AOGC-NGS RSGB_AML TGCM-AML 
##      323       45       90 
    
##     AMAS AOGC-NGS     MODY RSGB_AML     SDDS     SKDP TGCM-AML 
##        0       96       19       45       77       22       9

## pheno.ori[pheno.ori[, "projects"]=="AMAS",28:37]

 ## table(pheno.ori[, "Project"])

##   0   1 
## 220 135 
## sum(pheno[,"PD"]) ## pehno used to do asssocation tests - excluded PD
## sum(pheno.ori[,"PD"])
## sum(pheno[,"SCC"])





############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE
############################## ADD WEIGHT MODEL HERE


######################################################################################################
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##############    imafs begins    ####################################################################
######################################################################################################
#AN: commented out : for (ido.mafs in 1:length(do.MAFS)){


snpinfo.sliding<-FALSE
snpinfo.ori.keep<-snpinfo.ori ### cause sliding window replaces snpinfo.ori<-snpinfo.ori.keep
#do.MAFS<-c(0.001) # do.MAFS<-c(0.001,0.01) # do.MAFS<-c(0.001,0.005,0.01)
do.MAFS<-c(0.001)
ido.mafs<-1
#for (ido.mafs in 1:length(do.MAFS)){
    
######################################################################################################
n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.Control"]))
#n<-max(as.integer(summary.geno.extra[,"TOTAL.Alleles.AOGC"]))

    
p<-do.MAFS[ido.mafs]
# p<-0.01 ########### set MAF threshols HEREX1
# p<-0.005 ########### set MAF threshols HEREX1
# p<-0.05 ########### set MAF threshols HEREX1
    
sd.thresh<-4
n
p

alt.counts.thresh<-1
while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
alt.counts.thresh

## alt.counts.thresh<-1
## while( (alt.counts.thresh- n*p) / sqrt(n*p*(1-p)) <= sd.thresh){alt.counts.thresh<-alt.counts.thresh+1}
## alt.counts.thresh

n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.Control"])
# n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.AOGC"])
#  n  
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))

## alt.counts.thresh[1:50]
## summary.geno.extra[1:5,]


rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])<= alt.counts.thresh
rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh

## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC"])<= alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC.filt"])<= alt.counts.thresh
    
alt.counts.thresh.4.rare.in.Controls<-alt.counts.thresh

    
sum(rare.in.Control)
sum(rare.in.Control.filt )
names(rare.in.Control)<-key
names(rare.in.Control.filt )<-key


################################################ AOGC controls ################
if(sum( c("TOTAL.Alleles.AOGC","ALT.Alleles.AOGC","ALT.Alleles.AOGC.filt") %in%  colnames(summary.geno.extra))==3  ){
n<-as.numeric(summary.geno.extra[,"TOTAL.Alleles.AOGC"])
#  n  
alt.counts.thresh  <- round((sd.thresh*sqrt(n*p*(1-p)) + n*p))

alt.counts.thresh[1:50]
summary.geno.extra[1:5,]


## rare.in.Control<-as.numeric(summary.geno.extra[,"ALT.Alleles.Control"])<= alt.counts.thresh
## rare.in.Control.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.Control.filt"])<= alt.counts.thresh

rare.in.AOGC<-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC"])<= alt.counts.thresh
rare.in.AOGC.filt <-as.numeric(summary.geno.extra[,"ALT.Alleles.AOGC.filt"])<= alt.counts.thresh    
#alt.counts.thresh.4.rare.in.Controls<-alt.counts.thresh

    
sum(rare.in.AOGC)
sum(rare.in.AOGC.filt )
names(rare.in.AOGC)<-key
names(rare.in.AOGC.filt )<-key
}else{
rare.in.AOGC<-rep(TRUE,times=dim(a.indel)[1])
rare.in.AOGC.filt <- rep(TRUE,times=dim(a.indel)[1])
names(rare.in.AOGC)<-key
names(rare.in.AOGC.filt )<-key
}
################################################################################


maf.lt.all[1:5,]
maf.lt.all[1:5,]
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.001"])
#maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.01"])
maf.col<-paste("MAF.lt",p,sep=":")
maf.col
 maf.filter<-as.logical(maf.lt.all[,maf.col])
## maf.filter<-as.logical(maf.lt.all[,"MAF.lt:0.5"])
sum(maf.filter)
names(maf.filter)<-key


#pass<- rare.in.group & !no.genotypes & !high.missing & common.loci


## test<-"chr5:68667282:68667282:T:G:snp:68667282"
## summary.geno.extra[test,"ALT.Alleles.Control"]
## rare.in.Control[test]
######################################################################################################
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##############    itypes begins    ####################################################################
######################################################################################################


types<-c("Single Point","coding","non-coding","bad.effect","sliding.window")

#types<-c("coding")

itypes<-2



       
for(itypes in 1:length(types)){
#setwd("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/output_Analyse_project_MERGE_LEO_to_PIT_OneChr_showed_by_paul_before_genotype_reovery_line/")
  
snap.file<-types[itypes]
snap.file
#snap.file<-paste(snap.file,"QC_no27.FINAL_6sd",sep="")
## snap.file<-paste(snap.file,".withGQ.filter.recovery_weights_FULLQC_TIGHT_sd6_rare",sep="")
## snap.file<-paste(snap.file,".withGQ.filter.recovery_weights_TIGHT",sep="")
## snap.file<-paste(snap.file,".AOGC.","ALL",".FINAL.PCA",sep="") snap.file<-paste("Achal.",snap.file,".",p,".",project.files[ichr],sep="")
# snap.file<-paste("NEWwBIASwQG.WITH.RECOVERY",snap.file,".",p,".",project.files[ichr],sep="")
snap.file<-paste(snap.file.prefix,snap.file,".",p,".",project.files[ichr],sep="")
snap.file
## paste(snap.file,".withGQ.filter",sep="")

## paste(snap.file,".noGQ.filter",sep="")
## a.indel.use<-a.indel
## summary.geno.extra.use<-summary.geno.extra

## a.indel<-a.indel.ori
## summary.geno.extra<-summary.geno.extra.ori

##  a.indel<-a.indel.use
## summary.geno.extra<-summary.geno.extra.use

#& !novoalign.PD.bias & !bwa.PD.bias & !trueSeq.PD.bias & !nextera.PD.bias
#& !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias


# test.nextera[figure],test.trueSeq[figure],test.bwa[figure],test.novalign[figure],test.PD.nextera[figure],test.PD.trueSeq[figure],test.PD.bwa[figure],test.PD.novalign[figure],


if(types[itypes]=="Single Point"){
  
pass<- full.qual  & maf.filter & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  & GQ.Control.pass & GQ.cancer.pass

pass.all.cohorts<-    full.qual & maf.filter & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !in.common.hit.gene  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  & ok.missing & !very.high.missing.controls  & rare.in.Control.filt  & !are.repeats & !are.in.repeats &  GQ.Control.pass & GQ.cancer.pass

}

if(types[itypes]=="coding"){
  

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  & ok.missing & !very.high.missing.controls  & rare.in.Control.filt  & !are.repeats & !are.in.repeats & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias ## ORIGINALwBIAS

## pass.all.cohorts<- no.genotypes.cancer & !no.genotypes.PD & !no.genotypes.PD.filt & full.qual &  bad.coding  & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model  & ok.missing & !very.high.missing.controls  & rare.in.Control.filt  & !are.repeats & !are.in.repeats & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  ## ORIGINALwBIAS

    ## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias # NewwGias

    ## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  & GQ.Control.pass & GQ.cancer.pass # NEWwBIASwQG
    
    
pass<- full.qual &  bad.coding & maf.filter & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  & GQ.Control.pass & GQ.cancer.pass # NEWwBIASwQG & ok.missing  ## ver 9 -using & pass.possonian.control.model  no FALSE in thi set sum(!pass.possonian.control.model) & rare.in.AOGC

## pass<- full.qual &  bad.coding & maf.filter   & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt    & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias 
                                        #  & ok.missing  ## ver 9 -using & pass.possonian.control.model  no FALSE in thi set sum(!pass.possonian.control.model) & rare.in.AOGC

## pass.all.cohorts<- ( (!no.genotypes.PD & !no.genotypes.PD.filt) | (!no.genotypes.cancer ) ) &  full.qual &  bad.coding & maf.filter   & !in.common.hit.gene  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model  & !novoalign.bias & !bwa.bias & !trueSeq.bias & !nextera.bias  & GQ.Control.pass & GQ.cancer.pass


pass.all.cohorts<- full.qual  & !in.common.hit.gene & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  # ok.missing  ## ver 9 -using & pass.possonian.control.model  no FALSE in thi set sum(!pass.possonian.control.model) & rare.in.AOGC

sum(pass.all.cohorts)
#pass.all.cohorts<- full.qual &  bad.coding  & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model #
## sum(pass & !rare.in.AOGC.filt) 
## sum(pass & !rare.in.Control.filt)
}

if(types[itypes]=="bad.effect"){
  

pass<- full.qual &  bad.effect & maf.filter  & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt  & GQ.Control.pass & GQ.cancer.pass

pass.all.cohorts<- no.genotypes.cancer & !no.genotypes.PD & !no.genotypes.PD.filt & full.qual &  bad.effect  & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model  & ok.missing & !very.high.missing.controls  & rare.in.Control.filt  & !are.repeats & !are.in.repeats  & GQ.Control.pass & GQ.cancer.pass  #& ok.missing

}

if(types[itypes]=="non-coding"){
  

pass<- full.qual  & maf.filter & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & GQ.Control.pass & GQ.cancer.pass
pass.all.cohorts<- no.genotypes.cancer & !no.genotypes.PD & !no.genotypes.PD.filt & full.qual & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & bad.non.coding  & !in.common.hit.gene   & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model  & ok.missing & !very.high.missing.controls  & rare.in.Control.filt  & !are.repeats & !are.in.repeats & GQ.Control.pass & GQ.cancer.pass  #& ok.missing

}

if(types[itypes]=="sliding.window"){

pass<- full.qual  & maf.filter & !in.any.normal &  !in.any.normal.filt & !fail.bias.cancer & !fail.bias.control  & !in.common.hit.gene & !no.genotypes  & !unannotated.hits & not.flat.genotype  & hw.controls.ok.filt  &  rare.in.Control  & pass.possonian.aligner.model & pass.possonian.lib.model & !no.genotypes.filt & GQ.Control.pass & GQ.cancer.pass

sum(pass)

snpinfo.sliding<-sliding.window.generation(a.indel[pass,])
snpinfo.sliding[1:5,]

using.sliding.window<-TRUE
# snpinfo.ori.keep<-snpinfo.ori
snpinfo.ori<-snpinfo.sliding


}


##### full dump #  pass<-rep(TRUE,times=dim(a.indel)[1]) ; names(pass)<-key/media/TRI-T-DRIVE-tpleo/uqdi/UQCCG/Software/annovar/humandb/hg19/hg19_snp137_omim_maf.txt

sum(pass)
sum(pass)



## pass<-a.indel[,"Gene.Names"] %in% c("TP53","NOTCH1","NOTCH2","ATM","ACD","ASIP","BAP1","CASP8","CCND1","CDK4","MC1R","MITF","MTAP","MX2","OCA2","PARP1","PLA2G6","POT1","SLC45A2","TERF2IP","TERT","TYR","TYRP1","VDR","BCL2L12","KNSTRN","ISX","CDKN2A","BCL2L11","STK19","FJX1","TRHDE")
## sum(pass)


good.genotypes<-c(" ")
bad.genotypes<-c("chr2:240982245:240982245:-:GGT:indel")
#"chr12:49463361:49463376:TCCTCCACTCTTCCAT:-:indel" # RHEBL1 : rs200055056 Observed: -/CCTCCACTCTTCC includes AT where that last T delection is rare

# chr13:32753642:32753642:-:TGTGTGTGTGTG:indel:32753640 FRY dbSNP: rs200797182 Position: chr13:32753641-32753652 Observed: -/TGTGTGTGTGTG


# chr2:240982247:240982247:-:AGGGACGTGGGTGAAGAGCCGTGGGTGAAGGGCTGTGGGTGAAGAGCCGTGGG:indel this is the PRR21 indel locus that has two identical 
## chr2:240982245:240982245:-:GGT:indel this is the PRR21 indel locus that has two identical 
pass[bad.genotypes]


#AN
# chk.genotypes<-c("chr4:10099542:10099542:C:-:indel","chr4:10099540:10099540:C:-:indel")
# pass[chk.genotypes]
######
## a.indel[a.indel[,"start"]=="21854256",1:5]
## table[a.indel[,"Gene.Names"]
## chk<-a.indel[,"start"]=="21854256"
#AN
#summary.geno.extra[chk,]
#####

##       help[chk,]
## good.genotypes %in% names(pass)
## good.genotypes %in% snp.fail.filt
pass[ names(pass) %in% good.genotypes]<-TRUE
pass[ names(pass) %in% bad.genotypes]<-FALSE

    
help<-cbind( pass,full.qual,bad.coding,bad.non.coding,bad.effect,in.any.normal,in.any.normal.filt,fail.bias.cancer,fail.bias.control,in.common.hit.gene,unannotated.hits,not.flat.genotype,ok.missing,very.high.missing.controls,hw.controls.ok.filt,no.genotypes,no.genotypes.filt,is.benign.missense,pass.possonian.control.model,are.repeats,are.in.repeats,pass.possonian.aligner.model,pass.possonian.lib.model,bad.qual.locations, rare.in.Control , rare.in.Control.filt,GQ.Control.pass,GQ.cancer.pass,GQ.PD.pass)
#####################

## snpinfo.ori[1:5,]
## help[good.genotypes,]


##################### check for missing data



############################################
############################################
############################################
############### DO ASSOCIATION TEST  ##################
############################################
############################################


#############################

sum(pass)

genotypes<-a.indel[pass,the.samples.use] ## ordered correctly for phenotypes
snp.names<-key[pass] ## GEFOS ony name with start

#### snpinfo now A different size than a.indel since added pathways!!!  snpinfo[snpinfo[,"gene"]=="KNSTRN",]
snpinfo<-snpinfo.ori[snpinfo.ori[,"Name"] %in% snp.names,]

if(!exists("gene.weights")){
  gene.weights.subset<-1
}else{
gene.weights.subset<-gene.weights[snpinfo.ori[,"Name"] %in% snp.names] # weight in same order as snpinfo.ori
}
snpinfo<-cbind(snpinfo,gene.weights.subset)
#snpinfo[1:5,]
sum(is.na(as.numeric(snpinfo[,"gene.weights.subset"])))

###################################################/media/scratch/software/matlab/network.lic



if( sum(!(snp.names %in% snpinfo.ori[,"Name"]))>0){print("WARINING snp.names not in snpinfo- unusual!")}
dim(snpinfo)
length(snp.names)
dim(genotypes)


dim(genotypes)

print("start QC")
genotypes[genotypes=="NA"]<-NA
genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/1"]<-2

########### prevent any averaging
dim(genotypes)
genotypes[is.na(genotypes)]<-0
dim(genotypes)
########### prevent any averaging

########################################## GEFOS MINOR ALLELE TRANSFORMATION
## flip.geno<-gsub("2","3",genotypes[to.flip,])
## #flip.geno[1:15,1:10]
## flip.geno<-gsub("0","2",flip.geno)
## flip.geno<-gsub("3","0",flip.geno)
## genotypes[to.flip,]<-flip.geno
##########################################################################



#tapply(as.vector(genotypes),as.vector(genotypes),length)
num.col<-dim(genotypes)[2]
num.row<-dim(genotypes)[1]
## genotypes[1:5,1:20]
genotypes<-as.numeric(as.matrix(genotypes))
dim(genotypes)<-c(num.row,num.col)
genotypes<-t(genotypes) # samples x SNPS

colnames(genotypes)<-snp.names
rownames(genotypes)<-gsub(".GT$","",the.samples.use)



########## make snpInfo
#tapply(gene.list[,"CHR"],gene.list[,"CHR"],length)
## gene.list[gene.list[,"CHR"]=="Y",]
## a.indel[the.snps,1:20][1:12,]

### single pont tests: cohort.seq.single<-cohort.seq meta.results.burden.single<-meta.results.burden
### single pont tests:
## snpinfo.ori[1:5,]
## cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="gene",family=binomial(),verbose=FALSE)
## meta.results.burden<-singlesnpMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="gene")
### single pont tests:
### single pont tests:
### single pont tests:

if(types[itypes]=="Single Point"){
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="gene",family=binomial(),verbose=FALSE)
meta.results.burden<-singlesnpMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="gene")
meta.results.skat<-{}
meta.results.skatO<-meta.results.burden ### single point skat0 is meaningless so juat make the burden test
}else{

if(target.pheno.col %in% case.control){
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=binomial(),verbose=FALSE)
}else{
cohort.seq <- skatCohort(Z=genotypes,formula, SNPInfo = snpinfo, data=pheno,aggregateBy="cluster",family=gaussian(),verbose=FALSE) ## genes and clusters
}

meta.results.burden<-burdenMeta(cohort.seq,wts="gene.weights.subset",mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")

# meta.results.burden<-burdenMeta(cohort.seq,wts=1,mafRange = c(0,1),SNPInfo = snpinfo,aggregateBy="cluster")
#meta.results.skat<-skatMeta(cohort.seq,SNPInfo = snpinfo,aggregateBy="cluster") 

meta.results.skatO<-skatOMeta(cohort.seq,burden.wts =1,SNPInfo = snpinfo,aggregateBy="cluster",method = "integration") ## use this
# meta.results.skatO<-meta.results.burden
}


the.order<-     order(meta.results.burden[,"p"])
sum(is.na(meta.results.burden[,"p"])) ## bad p-values shoudl not happen
meta.results.burden<-  meta.results.burden[the.order,]
meta.results.burden[1:50,]
meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,] #meta.results.burden.3sd<-meta.results.burden

meta.results.skat<-{}
meta.results.skatO<-{}#  /media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/pac_bio.read.significance_NEW.r
## ## the.order<-     order(meta.results.skat[,"p"])
## ## meta.results.skat<-  meta.results.skat[the.order,]
## ## meta.results.skat[1:50,]

## the.order<-     order(meta.results.skatO[,"p"])
## sum(is.na(meta.results.skatO[,"p"])) ## bad p-values shoudl not happen
## meta.results.skatO<-  meta.results.skatO[the.order,]
#meta.results.skatO[1:20,]

meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
## ## meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,]
## clusters[1:5,]

## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,1],]
## meta.results.burden[meta.results.burden[,"gene"] %in% clusters[,2],]


## meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,]
if(types[itypes]=="sliding.window"){
  posns<-match(meta.results.burden[,"gene"],snpinfo[,"cluster"])
  missing<-is.na(posns)
  the.gene<-snpinfo[posns,"gene"]
  meta.results.burden<-cbind(the.gene,meta.results.burden)

}
                           

setwd(analysis.dir) #meta.results.skat<-{}
#output.dir<- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/skat_test/output"
#output.dir<- "/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul"
#setwd(output.dir)
## getwd()
## bad.non.coding
## snap.file<-"coding.0.01.all.geno.all.filters"
## snap.file<-"coding.0.001.all.geno.all.filters_no.imput"

##writing

# print(paste("Burden","ALL",snap.file,"txt",sep="."))
#  write.table(meta.results.burden,file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.skatO,file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/SKATO","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

#  print(paste("Burden","ALL",snap.file,"txt",sep="."))
#  write.table(meta.results.burden,file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/HPC/PBSHOME/aneupane/2015-11-11-AML/aml_burden_test_output_novoaligned2/Burden","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.skatO,file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/HPC/PBSHOME/aneupane/2015-11-11-AML/aml_burden_test_output_novoaligned2/SKATO","ALL",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/HPC/PBSHOME/aneupane/2015-11-11-AML/aml_burden_test_output_novoaligned2/Burden","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
#  write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste("/media/TRI-T-DRIVE-taneupan/uqdi/HPC/PBSHOME/aneupane/2015-11-11-AML/aml_burden_test_output_novoaligned2/SkatO","clusters",snap.file,"txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

print(paste("Burden","ALL",snap.file,"txt",sep="."))
write.table(meta.results.burden,file=paste0(path.save,type.project,"_Burden_ALL_",snap.file,".txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO,file=paste0(path.save,type.project,"_SKATO_ALL",snap.file,".txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.burden[meta.results.burden[,"gene"] %in% clusters.wanted,],file=paste0(path.save,type.project, "_Burden_clusters_",snap.file,".txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
write.table(meta.results.skatO[meta.results.skatO[,"gene"] %in% clusters.wanted,],file=paste0(path.save,type.project,"_SKATO_clusters_",snap.file,".txt"),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)


annotations<-a.indel[,c(1:6,16,28,7,30,34,37:42,43,14,32,33)]




getwd()
#save(list=c("synonymous","no.genotypes.cancer","cellularity","cancer","PDs","Control.alt.counts", "normal.alt.counts","the.samples.use","gene.weights","gene.weights.subset","filt","snp.fail.filt","use.wieght","weights","core.ann","case.control","snpinfo.ori","formula","clusters","pheno.types","ipheno","clusters.wanted","p","meta.results.skat","meta.results.skatO","meta.results.burden","pheno","pheno.ori","target.pheno.col","snpinfo","pass","high.missing.table","a.indel","a.indel.ori","help","key","summary.geno.extra","summary.geno.extra.ori","full.qual","bad.coding","bad.effect","maf.filter","in.common.hit.gene","on.x.y","unannotated.hits","not.flat.genotype","are.repeats","are.in.repeats","ok.missing","hw.controls.ok.filt","no.genotypes","rare.in.Control","rare.in.Control.filt","in.any.normal","in.any.normal.filt","are.in.repeats.back","are.in.repeats.forward","all.genes","contaminated","is.benign.missense","pass.possonian.control.model","pass.possonian.aligner.model","bad.qual.locations","filt"),file=paste(snap.file,".small_final.RData",sep="") )
print(paste("Got BURDEN/SKAT0 for ",project.files[ichr], ". Do the untangle",sep=""))

#source("/media/scratch/UQCCG_Pipeline_Rscripts/genotype.components.4.skatmeta.r") #/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/genotype.components.4.skatmeta.r")
source("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/genotype.components.4.skatmeta.r")

if(types[itypes]=="sliding.window"){
  snpinfo.ori<-snpinfo.ori.keep
  using.sliding.window<-FALSE
}



} # itype

# source("/media/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/genotype.components.4.skatmeta.r")

print(paste("Done: ",pheno.types[ipheno],"->",project.files[ichr]))
#save.image(file=paste(paste0("/media/TRI-T-DRIVE-taneupan/uqdi/UQCCG/UQCCG-Projects/Achal/Genetic_QC_findings/2015-08-14-AML/burden_test_code_new_by_paul/",snap.file,".RData",sep="")))
#save.image(file=paste(paste0("/media/TRI-T-DRIVE-taneupan/uqdi/HPC/PBSHOME/aneupane/2015-11-11-AML/aml_burden_test_output_novoaligned2/",snap.file,".RData",sep="")))
save.image(file=paste0(path.save,snap.file,"_",type.project,".RData"))

} # loop over imaf

  } # loop over projects ichr
  
  
} # loop over fam

