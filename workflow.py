import wolf
import pandas as pd
import dalmatian as dm

from wolF import *

def catchthefish(id,
                 gs_clean_bam,
                 gs_clean_bai,
                 OncoBed="gs://jba-utils/GRCh38-oncogenes-regions-CoMMPASS-adjusted-BCL.bed",
                 IgBed="gs://jba-utils/GRCh38-ig-regions-CoMMPASS-adjusted.bed",
                 min_map_quality_onco=60,
                 min_map_quality_ig=0,
                 exclude_flag_split=1540,
                 exclude_flag_mate=1548,
                 include_flag_mate=1,
                 genome="hg38",
                 eps=500,
                 minPts=2,
                 add_gs_suffix=True):
    
    bam_in_the_cloud = extract_sam(inputs={
                                        "id":id,
                                        "gs_clean_bam":gs_clean_bam,
                                        "gs_clean_bai":gs_clean_bai,
                                        "OncoBed":OncoBed,
                                        "IgBed":IgBed,
                                        "min_map_quality_onco":min_map_quality_onco,
                                        "min_map_quality_ig":min_map_quality_ig,
                                        "exclude_flag_split":exclude_flag_split,
                                        "exclude_flag_mate":exclude_flag_mate,
                                        "include_flag_mate":include_flag_mate
                                      })
    
    tx = cluster_reads(inputs = {
        "id":id,
        "OncoBed":OncoBed,
        "IgBed":IgBed,
        "reads":bam_in_the_cloud["pseudo_sam"],
        "genome":genome,
        "eps":eps,
        "minPts":minPts
    })


##############
#############
##############

with wolf.Workflow(
  workflow = catchthefish,
  conf = { "clust_frac": 0.7 }, # if you want to use more machines in the elastic cluster
  common_task_opts = { "retry" : 3} # will retry every task up to 5 times
) as w:
    w.run(id="Ultra_RT_1",
          gs_clean_bam="gs://7-genome-sphere/hg38-bam/Ultra_RT_1.bam",
          gs_clean_bai="gs://7-genome-sphere/hg38-bam/Ultra_RT_1.bam.bai")
