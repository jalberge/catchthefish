# variables
trx.annot <- structure(
  list(
    IG_CHR = c(
      "chr1",
      "chr11",
      "chr12",
      "chr14",
      "chr16",
      "chr17",
      "chr2",
      "chr20",
      "chr22",
      "chr3",
      "chr4",
      "chr4",
      "chr6",
      "chr6",
      "chr6",
      "chr8",
      "chr8"
    ),
    IG_START = c(
      118219025L,
      68896086L,
      4090496L,
      105294903L,
      78569424L,
      43155800L,
      89036865L,
      38347286L,
      22916326L,
      169638036L,
      1804726L,
      153354955L,
      7727323L,
      41783994L,
      108884686L,
      128142719L,
      144397552L
    ),
    IG_END = c(
      118415891L,
      69996562L,
      4324140L,
      106839645L,
      79382306L,
      43471519L,
      89312768L,
      38762745L,
      23411876L,
      169854050L,
      1919728L,
      153619446L,
      8154971L,
      42003789L,
      109102859L,
      129539512L,
      144502601L
    ),
    name = c(
      "FAM46C",
      "CCND1",
      "CCND2",
      "IGH",
      "MAF",
      "MAP3K14",
      "IGK",
      "MAFB",
      "IGL",
      "SAMD7",
      "FGFR3 / NSD2",
      "FBXW7/TMEM154",
      "TXNDC5",
      "CCND3",
      "FOXO3",
      "MYC",
      "MAFA"
    )
  ),
  class = "data.frame",
  row.names = c(NA,-17L)
)

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



vaf_col <- colorRamp2(seq(0, 0.5, length.out = 9), colors = brewer.pal(9, "PuRd"), transparency = 0, space = "LAB")

pal <- c( brewer.pal(8, "Dark2"),  brewer.pal(8, "Set1"),  brewer.pal(8, "Set2"))

cluster.all.tx.from.chr <- function(chr,
                                    delly,
                                    eps=2E5,
                                    minPts=5, 
                                    min_noise_partner=1000,
                                    min_ig_partner=1000, 
                                    chromosomes=chromosomes, 
                                    min.unique.patients=5,
                                    min.median.MAPQ=31,
                                    min.median.PE=1,
                                    min.median.SR=1){
  tx <- extract.tx.from.delly(delly, ig_chr=chr, ig_start=1, ig_end=300E6, chromosomes=chromosomes)
  
  clusters <- cluster.tx(tx, eps=eps, minPts=minPts, min_noise_partner=min_noise_partner)
  
  clusters$N_Unique_Patients <- unlist(lapply(str_split( clusters$patients, pattern = ", "), function(x){ length(unique(str_extract(string = x, pattern = "MMRF_[0-9]{4}"))) }))
  
  # keep only one of both contigs
  tx <- tx %>% filter(IG_CHR<PARTNER_CHR)
  clusters <- clusters %>% 
    filter(IG_CHR < PARTNER_CHR) %>%
    filter(across(any_of("N_Unique_Patients"), ~.x >= min.unique.patients)) %>%
    filter(mMAPQ >= min.median.MAPQ & mPE >= min.median.PE & mSR >= min.median.SR)
  # thanks to https://stackoverflow.com/questions/45146688/execute-dplyr-operation-only-if-column-exists
  
  bed1 <- tx.to.bed1(clusters)
  bed2 <- tx.to.bed2(clusters)
  
  list(tx, clusters, bed1, bed2)
}
