# CHM13 files for AnnotSV.

I tried to generate the chm13 files for AnnotSV, and here are the files and code.

The most important (files must be downloaded in T2T-CHM13):

RefSeq gene annotations:
hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
chm13: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/annovar/hs1_curGenev5.20.txt.gz
chm13: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/annovar/hs1_refGene.txt.gz

ENSEMBL gene annotations:
hg38: http://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.chr.gtf.gz
chm13: https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz

GC content annotations
hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
chm13: 
I couldn't find the corresponding file, but I generated it with the following code.
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz; mkdir chroms; gunzip -dc chm13v2.0_maskedY_rCRS.fa.gz | awk '{if($0~/^>/){chr=$1;gsub(/>| */,"",chr)};{print $0 >"chroms/"chr".fa"}}';tar -zcvf chm13.chromFa.tar.gz chroms/

Repeated sequences annotations
Ok, available at UCSC

Segmental duplication annotations
Ok, available at UCSC

GAP annotations
Ok, available at UCSC

Cytoband
hg38: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
chm13: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed

Others (liftover from GRCh38 to T2T-CHM13):

ClinVar pathogenic SV annotations
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240215.vcf.gz
Known pathogenic SNV/indel annotations
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240215.vcf.gz;
ClinVar benign SV annotations
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240215.vcf.gz

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20241223.vcf.gz;
gunzip -dc clinvar_20241230.vcf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($1~/^#/){print $0}else if($1=="MT"){$1="chrM";print $0}else{print "chr"$0}}' | bgzip >clinvar_20241230.reformat.vcf.gz;
java -jar /opt/software/broadinstitute/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar LiftoverVcf --RECOVER_SWAPPED_REF_ALT true -CHAIN /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain -R chm13v2.0_maskedY_rCRS.fa -I clinvar_20241230.reformat.vcf.gz -O clinvar_20241230.reformat.lifted.nochrM.chm13.vcf.gz -REJECT clinvar_20241230.reformat.rejected.vcf.gz;
(gunzip -dc clinvar_20241230.reformat.lifted.nochrM.chm13.vcf.gz;gunzip -dc clinvar_20241230.reformat.rejected.vcf.gz | grep "^chrM" ) | bgzip >clinvar_20241230.reformat.lifted.chm13.vcf.gz;

Dosage sensitive genes/regions annotation (ClinGen)
ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
ftp://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv
Not dosage sensitive genes/regions annotation (ClinGen)
https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv
wget https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz;
wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv;
wget https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv;
gunzip -dc Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{ if($3~/gene/ && $9~/gene_name/){g=$9;gsub(/.*gene_name "/,"",g);gsub(/".*/,"",g);print g"\tchr"$1":"$4"-"$5} }' >Homo_sapiens-GCA_009914755.4-2022_07-genes.gene.bed;
awk 'BEGIN{FS="\t";OFS="\t"}{ if(FNR==NR){a[$1]=$2}else if($0~/^#/){print $0}else if(a[$1]~/chr/){$4=a[$1];print $0}else{print $0 >"ClinGen_gene_curation_list_GRCh38.failed.tsv"} }' Homo_sapiens-GCA_009914755.4-2022_07-genes.gene.bed ClinGen_gene_curation_list_GRCh38.tsv | sed 's/GRCh38 (hg38): GCF_000001405.36/chm13 (hs1): GCA_009914755.4/g' >ClinGen_gene_curation_list_chm13.lifted.tsv;
awk 'BEGIN{FS="\t";OFS="\t"}{ if($1!~/^#/ && $4!~/tbd/){b=$4;gsub(/[:-]/,"\t",b);print b"\t"$0} }' ClinGen_region_curation_list_GRCh38.tsv >ClinGen_region_curation_list_GRCh38.reformat.bed;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab ClinGen_region_curation_list_GRCh38.reformat.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain ClinGen_region_curation_list_GRCh38.reformat.lifted.chm13.bed ClinGen_region_curation_list_GRCh38.reformat.failed.bed;
(head -n 6 ClinGen_region_curation_list_GRCh38.tsv | sed 's/GRCh38 (hg38): GCF_000001405.36/chm13 (hs1): GCA_009914755.4/g' ;cat ClinGen_region_curation_list_GRCh38.reformat.lifted.chm13.bed | awk 'BEGIN{FS="\t";OFS="\t"}{$7=$1":"$2"-"$3;print $0}' | cut -f 4-26 ) >ClinGen_region_curation_list_chm13.lifted.tsv;

dbVarNR pathogenic SV annotations
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.pathogenic.tsv.gz;
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.pathogenic.tsv.gz;
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/insertions/GRCh38.nr_insertions.pathogenic.tsv.gz;
gunzip GRCh38.nr_deletions.pathogenic.tsv.gz;
cat GRCh38.nr_deletions.pathogenic.tsv | grep -v "^#" | awk '{print "chr"$0}' >GRCh38.nr_deletions.pathogenic.reformat.bed;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_deletions.pathogenic.reformat.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_deletions.pathogenic.reformat.lifted.chm13.bed GRCh38.nr_deletions.pathogenic.reformat.failed.bed;
(grep ^# GRCh38.nr_deletions.pathogenic.tsv | sed 's/GRCh38/chm13/g';cat GRCh38.nr_deletions.pathogenic.reformat.lifted.chm13.bed;grep ^mt GRCh38.nr_deletions.pathogenic.tsv | sed 's/^mt/chrM/g') >chm13.nr_deletions.pathogenic.lifted.tsv;
gunzip GRCh38.nr_duplications.pathogenic.tsv.gz;
cat GRCh38.nr_duplications.pathogenic.tsv | grep -v "^#" | awk '{print "chr"$0}' >GRCh38.nr_duplications.pathogenic.reformat.bed;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_duplications.pathogenic.reformat.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_duplications.pathogenic.reformat.lifted.chm13.bed GRCh38.nr_duplications.pathogenic.reformat.failed.bed;
(grep ^# GRCh38.nr_duplications.pathogenic.tsv | sed 's/GRCh38/chm13/g';cat GRCh38.nr_duplications.pathogenic.reformat.lifted.chm13.bed;grep ^mt GRCh38.nr_duplications.pathogenic.tsv | sed 's/^mt/chrM/g') >chm13.nr_duplications.pathogenic.lifted.tsv;
gunzip GRCh38.nr_insertions.pathogenic.tsv.gz;
cat GRCh38.nr_insertions.pathogenic.tsv | grep -v "^#" | awk 'BEGIN{FS="\t";OFS="\t"}{$3=$3+1;print "chr"$0}' >GRCh38.nr_insertions.pathogenic.reformat.bed;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_insertions.pathogenic.reformat.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_insertions.pathogenic.reformat.lifted.chm13.bed GRCh38.nr_insertions.pathogenic.reformat.failed.bed;
(grep ^# GRCh38.nr_insertions.pathogenic.tsv | sed 's/GRCh38/chm13/g';cat GRCh38.nr_insertions.pathogenic.reformat.lifted.chm13.bed |awk 'BEGIN{FS="\t";OFS="\t"}{$3=$3-1;print $0}';grep ^mt GRCh38.nr_insertions.pathogenic.tsv | sed 's/^mt/chrM/g') >chm13.nr_insertions.pathogenic.lifted.tsv;

gnomAD benign SV annotations
https://gnomad.broadinstitute.org/downloads#v4-structural-variants
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz
gunzip gnomad.v4.1.sv.sites.bed.gz
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab gnomad.v4.1.sv.sites.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain gnomad.v4.1.sv.sites.lifted.col1to64.noheader.chm13.bed  gnomad.v4.1.sv.sites.failed.bed
awk 'BEGIN{FS="\t";OFS="\t"}{ if(FNR==NR){a[$4]=$1"\t"$2"\t"$3}else{ if($1~/^#/){print $0}else if(a[$4]~/chr/){ split(a[$4],b,"\t");$1=b[1];$2=b[2];$3=b[3];print $0} } }' gnomad.v4.1.sv.sites.lifted.col1to64.noheader.chm13.bed gnomad.v4.1.sv.sites.bed | bgzip >gnomad.v4.1.sv.sites.lifted.col1to626.chm13.bed.gz

DGV benign SV annotations
http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt
cat GRCh38_hg38_variants_2020-02-25.txt | awk 'BEGIN{FS="\t";OFS="\t"}{ if($6~/insertion/ && $3==$4){$4=$4+1};if( $1!="variantaccession" && $2!="N" && $2!="" ){print "chr"$2,$3,$4,$0} }' >GRCh38_hg38_variants_2020-02-25.txt.reformat.bed
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38_hg38_variants_2020-02-25.txt.reformat.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38_hg38_variants_2020-02-25.txt.reformat.lifted.chm13.bed GRCh38_hg38_variants_2020-02-25.txt.reformat.failed.bed
(head -n 1 GRCh38_hg38_variants_2020-02-25.txt;cat GRCh38_hg38_variants_2020-02-25.txt.reformat.lifted.chm13.bed | awk 'BEGIN{FS="\t";OFS="\t"}{$5=$1;$6=$2;$7=$3;print $0}' | cut -f 4-23 ) >chm13_hs1_variants_2020-02-25.txt

DDD benign SV annotations
https://www.deciphergenomics.org/files/downloads/population_cnv_grch38.txt.gz
wget https://www.deciphergenomics.org/files/downloads/population_cnv_grch38.txt.gz
cat population_cnv_grch38.txt | awk 'BEGIN{FS="\t";OFS="\t"}{ if($1!~/^#/){print "chr"$2"\t"$3"\t"$4"\t"$0} }' >population_cnv_grch38.reformat.bed
(head -n 1 population_cnv_grch38.txt;cat population_cnv_grch38.reformat.lifted.chm13.bed | awk 'BEGIN{FS="\t";OFS="\t"}{$5=$1;$6=$2;$7=$3;print $0}' | cut -f 4-19 ) >population_cnv_chm13.txt

1000 genomes benign SV annotations
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz
gunzip -dc ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz >ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf
cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf | grep -v ^# | awk 'BEGIN{FS="\t";OFS="\t"}{ end=$2+1;if($8~/;END=/){end=$8;gsub(/.*;END=/,"",end);gsub(/;.*/,"",end)};print "chr"$1"\t"$2"\t"end"\t"$3 }' >ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.bed
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.lifted.chm13.bed ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.failed.bed
(cat ALL.wgs.mergedSV.v8.20130502.svs.genotypes.chm13.vcf.header;awk 'BEGIN{FS="\t";OFS="\t"}{ if(FNR==NR){a[$4]=$0}else if(a[$3]~/^chr/){split(a[$3],b,"\t");$1=b[1];$2=b[2];gsub(/;END=[0-9]*;/,";END="b[3]";",$8);print $0} }' ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.lifted.chm13.bed ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf ) | bgzip >ALL.wgs.mergedSV.v8.20130502.svs.genotypes.chm13.vcf.gz

dbVar benign SV annotations
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.common.bed.gz;
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.common.bed.gz;
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/insertions/GRCh38.nr_insertions.common.bed.gz;
gunzip GRCh38.nr_deletions.common.bed.gz;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_deletions.common.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_deletions.common.lifted.chm13.bed GRCh38.nr_deletions.common.failed.bed;
cat GRCh38.nr_deletions.common.lifted.chm13.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_del"}' >chm13.nr_deletions.common.lifted.bed;
gunzip GRCh38.nr_deletions.common.bed.gz;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_duplications.common.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_duplications.common.lifted.chm13.bed GRCh38.nr_duplications.common.failed.bed;
cat GRCh38.nr_duplications.common.lifted.chm13.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_dup"}' >chm13.nr_duplications.common.lifted.bed;
gunzip GRCh38.nr_insertions.common.bed.gz;
/raid5/DATABASE/UCSC/liftOver -bedPlus=3 -tab GRCh38.nr_insertions.common.bed /raid5/DATABASE/T2T_analysis/grch38-chm13v2.chain GRCh38.nr_insertions.common.lifted.chm13.bed GRCh38.nr_insertions.common.failed.bed;
cat GRCh38.nr_insertions.common.lifted.chm13.bed | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_ins"}' >chm13.nr_insertions.common.lifted.bed;
