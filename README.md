# bedGraph2Cluster
MATLAB Function for the k-Means Clustering of ChIP-seq bedGraph Data

| MATLAB_VERSION  | RUNTIME_VERSION |
| ------------- | ------------- |
| R2018b | 9.5 |

```MATLAB
bedGraph2Cluster(bedGraphs_Peak, bedGraphs_Cluster, bedGraphs_Noncluster, bedGraphs_Control, Outdir, BED_Bin, k, QNorm, Workingdir)
```
<br/>

## Required arguments
**bedGraphs_Peak (string)**: comma-delimited list of bedGraph files to be included during peak calling

**bedGraphs_Cluster (string)**: comma-delimited list of bedGraph files to be included during *k*-means clustering

**bedGraphs_Noncluster (string)**: comma-delimited list of bedGraph files to be excluded during *k*-means clustering

**bedGraphs_Control (string)**: comma-delimited list of bedGraph files to be used as controls for peak calling

**Outdir (string)**: path to the output directory

**BED_Bin (string)**: path to the BED file used for binned bedGraph generation

**k (string)**: number of clusters during *k*-means clustering

**QNorm (string)**: whether to perform QNorm normalization ("true": QNorm, "false": CPM)

## Optional arguments
**Workingdir (string)**: path to the output directory

<br/>

# Examplary usage
## with license for MATLAB
```MATLAB
bedGraph2Cluster("bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph", "bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph", "bedgraph/E2F1.bedgraph,bedgraph/CTCF.shSCR.bedgraph,bedgraph/c-Jun.shSCR.bedgraph", "bedgraph/INPUT.WT.bedgraph,bedgraph/INPUT.dCDK.bedgraph", "test_output", "bed/hg19.200bp.bed", "8", "true", "../")
```
## without license for MATLAB
[![PyPI version](https://badge.fury.io/py/run_matlab.svg)](https://badge.fury.io/py/run_matlab)
```shell
pip3 insatll run_matlab
run_matlab install -v R2018b -r 9.5
git clone https://github.com/hanjunlee21/bedGraph2Cluster
run_matlab run -v R2018b -r 9.5 bedGraph2Cluster bedGraph2Cluster bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph bedgraph/E2F1.bedgraph,bedgraph/CTCF.shSCR.bedgraph,bedgraph/c-Jun.shSCR.bedgraph bedgraph/INPUT.WT.bedgraph,bedgraph/INPUT.dCDK.bedgraph test_output bed/hg19.200bp.bed 8 true $PWD
```
