import wolf

class extract_sam(wolf.Task):
    name = "extract_sam"
    inputs = {
        "id":None,
        "gs_clean_bam":None,
        "gs_clean_bai":None,
        "OncoBed":None,
        "IgBed":None,
        "min_map_quality_onco":60,
        "min_map_quality_ig":0,
        "exclude_flag_split":1540,
        "exclude_flag_mate":1548,
        "include_flag_mate":1,
        "add_gs_suffix":True
    }
    overrides = {"gs_clean_bam":"string", "gs_clean_bai":"string"}
    script = """
    set -euxo pipefail

    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

    # in case we need to reprocess a bit the URI
    bam="${gs_clean_bam}"
    bai="${gs_clean_bai}"

    cat ${OncoBed} | \
        awk '{ print "$7==\\""$1"\\"","&& $8>"$2,"&& $8<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_MP_Onco.awk

    cat ${OncoBed} | \
        awk '{ print "$12==\\""$1"\\"","&& $13>"$2,"&& $13<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_SA_Onco.awk

    cat ${IgBed} | \
        awk '{ print "$7==\\""$1"\\"","&& $8>"$2,"&& $8<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_MP_Ig.awk

    cat ${IgBed} | \
        awk '{ print "$12==\\""$1"\\"","&& $13>"$2,"&& $13<="$3}' | \
        sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n/ || /g' | \
        cat <(echo -n "{ if (") - | \
        awk '{print $0" ) { print $0 ;}}"}' > \
        code_filter_sam_SA_Ig.awk

    samtools view \
        -q ${min_map_quality_ig} \
        -F ${exclude_flag_mate} \
        -f ${include_flag_mate} \
        -M -L ${IgBed} -X ${bam} ${bai} | \
        cut -f1-11 | \
        awk -f code_filter_sam_MP_Onco.awk | \
        awk '{print "Paired-Read",$0}' OFS='\t' \
        > ${id}_MP_reads_Ig_Onco.txt

    samtools view \
        -q ${min_map_quality_ig} \
        -F ${exclude_flag_split} \
        -M -L ${IgBed} -e '[SA]' -X ${bam} ${bai} | \
        cut -f1-12 | \
        sed 's/SA:Z://g' | \
        awk -F'\t' 'gsub(/,/,"\t",$12)' OFS='\t' | \
        awk -f code_filter_sam_SA_Onco.awk | \
        cut -f1-11 | \
        awk '{print "Split-Read",$0}' OFS='\t' \
        > ${id}_SA_reads_Ig_Onco.txt

    samtools view \
        -q ${min_map_quality_onco} \
        -F ${exclude_flag_mate} \
        -f ${include_flag_mate} \
        -M -L ${OncoBed} -X ${bam} ${bai} | \
        cut -f1-11 | \
        awk -f code_filter_sam_MP_Ig.awk | \
        awk '{print "Paired-Read",$0}' OFS='\t' \
        > ${id}_MP_reads_Onco_Ig.txt

    samtools view \
        -q ${min_map_quality_onco} \
        -F ${exclude_flag_split} \
        -M -L ${OncoBed} -e '[SA]' -X ${bam} ${bai} | \
        cut -f1-12 | \
        sed 's/SA:Z://g' | \
        awk -F'\t' 'gsub(/,/,"\t",$12)' OFS='\t' | \
        awk -f code_filter_sam_SA_Ig.awk | \
        cut -f1-11 | \
        awk '{print "Split-Read",$0}' OFS='\t' \
        > ${id}_SA_reads_Onco_Ig.txt

    cat ${id}_MP_reads_Ig_Onco.txt ${id}_SA_reads_Ig_Onco.txt ${id}_MP_reads_Onco_Ig.txt ${id}_SA_reads_Onco_Ig.txt | \
    awk -v id=${id} -F'\t' '{print id,$0}' OFS='\t' > ${id}_tmp_reads.txt

    cut -f3 ${id}_tmp_reads.txt | sort | uniq > potential_reads.txt
    cat ${OncoBed} ${IgBed} > potential_regions.txt

    samtools view \
        -o "${id}_potential_reads.bam" \
        -F ${exclude_flag_split} \
        -M -L potential_regions.txt \
        -N potential_reads.txt \
        -X ${bam} ${bai}
    samtools index "${id}_potential_reads.bam"

    echo -e "ID\tREADTYPE\tQNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tQUAL" | cat - "${id}_tmp_reads.txt" > "${id}_mini_reads.txt"

    """
    outputs = {  "pseudo_sam" : "*_mini_reads.txt",
        "potential_reads_bam" : "*_potential_reads.bam",
        "potential_reads_bai" : "*_potential_reads.bam.bai"
              }
    docker = "jbalberge/samtools_cloud:1.13"


class cluster_reads(wolf.Task):
    name = "cluster_reads"
    inputs = {
        "id":None,
        "OncoBed":None,
        "IgBed":None,
        "reads":None,
        "genome":"hg19",
        "eps":500,
        "minPts":2
    }
    script = """
    set -euxo pipefail

    Rscript - \
    "${reads}" \
    "${OncoBed}" \
    "${IgBed}" \
    "${genome}" \
    "${eps}" \
    "${minPts}" \
    "${id}_clusters_eps${eps}_minPts${minPts}_raw.txt" \
    "${id}_clusters_eps${eps}_minPts${minPts}_filtered.txt" \
    "${id}_clusters_eps${eps}_minPts${minPts}.pdf" << "HEREDOC"

sem <- function(x) sd(x)/sqrt(length(x))

annotate.position <- function(hit.chr, hit.pos, bed) {
  name <- bed %>% filter(chr == hit.chr & start <= hit.pos & end >= hit.pos) %>% pull(name)
  if( identical(name, character(0))) {
    return("")
  } else {
    return(name)
  }
}


run.dbscan <- function(df, eps=1000, min.pts=2){
  mat <- df %>% select(POS.IG, POS.ONCO) %>% as.matrix()
  res <- dbscan(mat, eps, minPts = min.pts, weights = NULL, borderPoints = TRUE)
  df$Cluster=res$cluster
  df
}



extract.tx.from.delly <- function(delly.output,
                                  ig_chr="14",
                                  ig_start=105586000,
                                  ig_end=106880000,
                                  chromosomes=c(1:22, "X", "Y")) {
  delly.output %>%
    filter(SVTYPE=="TRA") %>%
    filter( ( CHROM==ig_chr & POS >= ig_start & POS <= ig_end ) | ( CHR2 == ig_chr & ENDPOSSV >= ig_start & ENDPOSSV <= ig_end )) %>%
    mutate(PARTNER_CHR = ifelse(CHROM==ig_chr, CHR2, CHROM),
           PARTNER_POS = ifelse(CHROM==ig_chr, ENDPOSSV, POS),
           IG_CHR = ig_chr,
           IG_POS = ifelse(CHROM==ig_chr, POS, ENDPOSSV)) %>%
    filter(CHR2 %in% chromosomes & CHROM %in% chromosomes) %>%
    rowwise() %>%
    mutate(VAF=if(PRECISE=="true") TUMOR_RV / ( TUMOR_RV + TUMOR_RR ) else ( TUMOR_DV / ( TUMOR_DV + TUMOR_DR )) )
}


cluster.tx <- function(tx, eps=2E5, minPts=3, min_noise_partner=1000, min_ig_partner=1000) {
  tx %>%
    group_by(IG_CHR, PARTNER_CHR) %>%
    group_modify(~ {
      distances <-  dist(matrix(c(.x$PARTNER_POS, .x$IG_POS), ncol=2, byrow = FALSE), method = "maximum");
      .x %>% mutate(Cluster = dbscan(distances, eps=eps, minPts=minPts)$cluster)}) %>%
    mutate(Chr_Cluster=paste0(IG_CHR, "-", PARTNER_CHR, "_", Cluster)) %>%
    filter(Cluster!=0) %>%
    group_by(Chr_Cluster, IG_CHR, PARTNER_CHR) %>%
    summarise(COUNT_CLUSTERED_TX=n(),
              IG_START=min(IG_POS),
              IG_END=max(IG_POS),
              PARTNER_START=min(PARTNER_POS),
              PARTNER_END=max(PARTNER_POS),
              INTERVAL_LENGTH=PARTNER_END-PARTNER_START,
              IG_LENGTH=IG_END-IG_START,
              patients=paste(Specimen_ID, collapse = ", "),
              mMAPQ=median(MAPQ),
              mPE=median(as.numeric(PE), na.rm = TRUE),
              mSR=median(as.numeric(SR), na.rm = TRUE),
              mVAF=median(VAF), na.rm = TRUE) %>%
    filter(INTERVAL_LENGTH > min_noise_partner & IG_LENGTH >= min_ig_partner) %>%
    ungroup()
}

tx.to.bed1 <- function(x) {
  x %>%
    mutate(IG_CHR=paste0("chr", IG_CHR)) %>%
    dplyr::select(IG_CHR, IG_START, IG_END)
}
tx.to.bed2 <- function(x) {
  x %>%
    mutate(PARTNER_CHR=paste0("chr", PARTNER_CHR)) %>%
    dplyr::select(PARTNER_CHR, PARTNER_START, PARTNER_END)
}

# This script takes potential reads from CatchTheFISH,
# pairs reads, and cluster them, to call translocations

# INIT --------------------------------------------------------------------

set.seed(123)

library(tidyverse)
library(dbscan)
library(circlize)

args = commandArgs(trailingOnly=TRUE)

# INPUT
reads <- args[1]
onco.bed <- args[2]
ig.bed <- args[3]

# PARAMS
genome <- as.character(args[4])
eps <- as.numeric(args[5])
minPts <- as.numeric(args[6])

# OUTPUT
clustered.tx <- args[7]
clustered.tx.filtered <- args[8]
circos.plot <- args[9]

# Parse input files -------------------------------------------------------

tx.header <- c("ID","SVTYPE","READNAME","FLAG","CHR1","POS1","QUAL","CIGAR","CHRB","POSB","TLEN","SEQ","SEQQUAL")
tx.col.classes <- c("character", "character", "character", "numeric", "character", "numeric", "numeric", "character", "character", "numeric", "numeric", "character", "character")
tx <- read.delim(reads, col.names = tx.header, colClasses = tx.col.classes, header = TRUE)

if( nrow(tx ) < 2 ) {
  # first, in case there is no abnormal read called / potential translocation
  res.colnames <- c("ID","Chr_Cluster","IG","ONCO","CHR.IG","IG_START","IG_END","CHR.ONCO","PARTNER_START","PARTNER_END","INTERVAL_LENGTH","IG_LENGTH","reads","maxMAPQ.ONCO","maxMAPQ.IG","N_Matches","N_Molecules","N_Mate_Pairs","N_Split_Reads","Breakpoint")
  res <- setNames(data.frame(matrix(ncol = length(res.colnames), nrow = 0)), res.colnames)
  pdf(circos.plot, paper = "a4")
  circos.clear()
  circos.par("track.height" = 0.1, start.degree = 90)
  circos.initializeWithIdeogram(species = genome,  track.height = convert_height(2, "mm"), c("ideogram", "labels"))
  dev.off()
} else {

  bed.header <- c("chr", "start", "end", "name")
  bed.col.classes <- c("character", "numeric", "numeric", "character")
  onco.bed <- read.delim(onco.bed, col.names = bed.header, colClasses = bed.col.classes, header = FALSE)
  ig.bed <- read.delim(ig.bed, col.names = bed.header, colClasses = bed.col.classes, header = FALSE)

  tx.2 <- tx %>%
    rowwise() %>%
    mutate(IG=annotate.position(hit.chr = CHR1, hit.pos = POS1, bed = ig.bed),
           ONCO=annotate.position(hit.chr = CHR1, hit.pos = POS1, bed = onco.bed)) %>%
    group_by(ID, SVTYPE, READNAME) %>%
    filter(n()>=2)
  # if read name appears only once, it could be that mate was annotated invalid by vendor

  tx.ig <- tx.2 %>%
    filter(IG!="") %>%
    transmute(ID, SVTYPE, READNAME, CHR.IG=CHR1, POS.IG=POS1, IG, QUAL.IG=QUAL)

  tx.onco <- tx.2 %>%
    filter(ONCO!="") %>%
    transmute(ID, SVTYPE, READNAME, CHR.ONCO=CHR1, POS.ONCO=POS1, ONCO, QUAL.ONCO=QUAL)

  # should be only one row per read / SVTYPE / patient
  # QUAL ONCO is more important than IG (we want to be very stringent on oncogenic partner, while IG is ~~~)
  tx.join <- inner_join(tx.ig, tx.onco, by=c("ID", "SVTYPE", "READNAME")) %>%
    group_by(ID, SVTYPE, READNAME) %>%
    slice_max(n = 1, order_by = QUAL.ONCO, with_ties = FALSE)

  tx.join.annot <- tx.join %>%
    left_join(onco.bed, by=c("ONCO"="name"), suffix=c("", "oncoBed")) %>%
    left_join(ig.bed, by=c("IG"="name"), suffix=c("", "_igBed")) %>%
    mutate(trx_name=
             paste0("t(",
                    min(as.numeric(CHR.IG), as.numeric(CHR.ONCO)),
                    ";",
                    max(as.numeric(CHR.IG), as.numeric(CHR.ONCO)),
                    ")"),
           trx_full_name=paste0(IG, "-", ONCO))

  # Clustering --------------------------------------------------------------

  res1 <- tx.join.annot %>%
    group_by(ID, ONCO, IG) %>%
    group_modify(~ {
      distances <-  dist(matrix(c(.x$POS.IG, .x$POS.ONCO), ncol=2, byrow = FALSE), method = "maximum");
      .x %>% mutate(Cluster = dbscan(distances, eps=eps, minPts=minPts)$cluster)}) %>%
    mutate(Chr_Cluster=paste0(IG, "-", ONCO, "_", Cluster)) %>%
    filter(Cluster!=0) %>%
    group_by(ID, Chr_Cluster, IG, ONCO) %>%
    summarise(CHR.IG=first(CHR.IG),
              IG_START=min(POS.IG),
              IG_END=max(POS.IG),
              CHR.ONCO=first(CHR.ONCO),
              PARTNER_START=min(POS.ONCO),
              PARTNER_END=max(POS.ONCO),
              INTERVAL_LENGTH=PARTNER_END-PARTNER_START,
              IG_LENGTH=IG_END-IG_START,
              reads=paste(READNAME, collapse = ","),
              MAPQ.ONCO=paste(QUAL.ONCO,  collapse = ","),
              MAPQ.IG=paste(QUAL.IG,  collapse = ","),
              maxMAPQ.ONCO=max(QUAL.ONCO),
              maxMAPQ.IG=max(QUAL.IG),
              N_Matches=n(),
              N_Molecules=length(unique(READNAME)),
              N_Mate_Pairs=sum(SVTYPE=="Paired-Read"),
              N_Split_Reads=sum(SVTYPE=="Split-Read"),
              Breakpoint=paste0("Breakpoint_", first(Cluster))) %>%
    ungroup()

  # Further filtering -------------------------------------------------------

  res <- res1 %>%
      filter(MAPQ.ONCO>=55 & INTERVAL_LENGTH>=1 & IG_LENGTH>=1)

  # Circos plot -------------------------------------------------------------

  pdf(circos.plot, paper = "a4")
  circos.clear()
  circos.par("track.height" = 0.1, start.degree = 90)
  circos.initializeWithIdeogram(species = genome,  track.height = convert_height(2, "mm"), c("ideogram", "labels"))

  if(nrow(res)>=1) {

    annot.partners <- res %>%
      transmute(chr=ifelse(startsWith(CHR.ONCO, "chr"), CHR.ONCO, paste0("chr", CHR.ONCO)),
                start=as.numeric(PARTNER_START),
                end=as.numeric(PARTNER_END),
                name=as.character(ONCO))

    annot.ig <- res %>%
      transmute(chr=ifelse(startsWith(CHR.IG, "chr"), CHR.IG, paste0("chr", CHR.IG)),
                start=as.numeric(IG_START),
                end=as.numeric(IG_END),
                name=as.character(IG))

    circos.genomicTrack(annot.partners, ylim=c(0, 0.1),
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ytop = 0.1, ybottom = 0,...)
                          circos.genomicText(region, value, y=0, labels.column = 1, facing = "clockwise", niceFacing = TRUE, cex=0.8, adj=degree(0), ...)
                        }, bg.border = NA )

    circos.genomicLink(region1 = annot.partners,
                       region2 = annot.ig,
                       col = "red",
                       lwd = 10)
  }

  title(paste(unique(res$ID), collapse=", "))
  dev.off()
}


# Write output ------------------------------------------------------------

write.table(res, clustered.tx.filtered, row.names = FALSE, quote = TRUE, sep = "\t")
write.table(res1, clustered.tx, row.names = FALSE, quote = TRUE, sep = "\t")

HEREDOC

    """
    outputs = {
        "clusters_tx" : "*_raw.txt",
        "clusters_tx_filtered" : "*_filtered.txt",
        "circos_tx" : "*.pdf"
    }
    docker = "jbalberge/r-dbscan:latest"


