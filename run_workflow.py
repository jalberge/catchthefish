import wolf
import pandas as pd
import dalmatian as dm

from wolF import *

def terra_na(x):
    y = ( pd.isna(x) ) | ( x=="")
    return y

WORKSPACE = "broad-getzlab-mm-germline-t/MM_WGS_GenomeSphere"

# find the corresponding pairs
WIC = wolf.fc.WorkspaceInputConnector(WORKSPACE)

P = WIC.get_pairs_as_joint_samples()

P = P.merge(WIC.pairs, left_index=True, right_index=True, suffixes=[None, '_pair'])

P = P.loc[ ~terra_na(P["hg38_analysis_ready_bam_T"]) & ( ( terra_na(P["clusters_tx"] ) if "clusters_tx" in P.columns else [ True for x in P.index ] ) ) ]

#############
##############

with wolf.Workflow(
  workflow = catchthefish,
  conf = { "clust_frac": 0.7 }, # if you want to use more machines in the elastic cluster
  common_task_opts = { "retry" : 3} # will retry every task up to 5 times
) as w:
    for pair,p in P.iterrows():
        w.run(id=pair,
                run_name=pair,
                gs_clean_bam=p["hg38_analysis_ready_bam_T"],
                gs_clean_bai=p["hg38_analysis_ready_bam_index_T"],
                workspace = WORKSPACE)
