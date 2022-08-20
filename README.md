# bedGraph2Cluster
MATLAB Function for the k-Means Clustering of ChIP-seq bedGraph Data

| MATLAB_VERSION  | RUNTIME_VERSION |
| ------------- | ------------- |
| R2018b | 9.5 |

```MATLAB
bedGraph2Cluster(bedGraphs_Target, bedGraphs_Nontarget, bedGraphs_Control, Outdir, BED_Bin, k, QNorm, Workingdir)
```
<br/>

## Required arguments
**bedGraphs_Target (string)**: comma-delimited list of bedGraph files to be included during *k*-means clustering

**bedGraphs_Nontarget (string)**: comma-delimited list of bedGraph files to be excluded during *k*-means clustering

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
bedGraph2Cluster("bam/RB.WT.filtered.bedgraph,bam/RB.dCDK.filtered.bedgraph", "bam/E2F1.filtered.bedgraph,bam/CTCF.shSCR.filtered.bedgraph,bam/c-Jun.shSCR.filtered.bedgraph", "bam/INPUT.WT.filtered.bedgraph,bam/INPUT.dCDK.filtered.bedgraph", "test_output", "bed/hg19.200bp.bed", "8", "true", "../")
```
## without license for MATLAB
[![PyPI version](https://badge.fury.io/py/run_matlab.svg)](https://badge.fury.io/py/run_matlab)
```shell
pip3 insatll run_matlab
run_matlab install -v R2018b -r 9.5
git clone https://github.com/hanjunlee21/bedGraph2Cluster
run_matlab run -v R2018b -r 9.5 bedGraph2Cluster bedGraph2Cluster bam/RB.WT.filtered.bedgraph,bam/RB.dCDK.filtered.bedgraph bam/E2F1.filtered.bedgraph,bam/CTCF.shSCR.filtered.bedgraph,bam/c-Jun.shSCR.filtered.bedgraph bam/INPUT.WT.filtered.bedgraph,bam/INPUT.dCDK.filtered.bedgraph test_output bed/hg19.200bp.bed 8 true $PWD
```
