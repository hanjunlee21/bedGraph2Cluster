# bedGraph2Cluster
**MATLAB Function for the k-Means Clustering of ChIP-seq bedGraph Data**

[![DOI](https://zenodo.org/badge/526668850.svg)](https://zenodo.org/badge/latestdoi/526668850)

<br/>

| MATLAB_VERSION  | RUNTIME_VERSION |
| ------------- | ------------- |
| R2018b | 9.5 |

```MATLAB
bedGraph2Cluster(bedGraphs_Signal, bedGraphs_Control, bedGraphs_Cluster, bedGraphs_Heatmap, outdir, bed_bin, fold_change, normalization_method, k, distance_method, clustering_method, workingdir)
```

<br/>

## Required arguments
**bedGraphs_Signal (string)**: comma-delimited list of bedGraph files to be included during peak calling

**bedGraphs_Control (string)**: comma-delimited list of bedGraph files to be used as controls for peak calling

**bedGraphs_Cluster (string)**: comma-delimited list of bedGraph files to be included during *k*-means clustering

**bedGraphs_Heatmap (string)**: comma-delimited list of bedGraph files to be included in heatmap

**outdir (string)**: path to the output directory

**bed_bin (string)**: path to the BED file used for binned bedGraph generation

**fold_change (string)**: threshold for the fold change over the control during peak calling

**normalization_method (string)**: normalization method to utilize ("QNorm": QNorm, "CPM": CPM)

**k (string)**: number of clusters during *k*-means clustering

**distance_method (string)**: distance metric for *k*-means clustering ("sqeuclidean", "cityblock", "cosine", "correlation")

**clustering_method (string)**: clustering method to utilize ("1": profile, "2": profile+scalar, "3": symmetry_collapsed_profile+scalar)

## Optional arguments
**workingdir (string)**: path to the output directory

## Output files
**tiles_200_data.mat**: MATLAB variable for the entire dataset

**tiles_200_data_peak.mat**: MATLAB variable for only the peakset

**peaks.bed**: BED file of the peakset

**peaks.clust${k}.${clusterID}.bed**: BED file of the peakset for each cluster

**clustering_heatmap.pdf**: heatmap for the *k*-means clustering

<br/>

# Examplary usage
## with license for MATLAB
```MATLAB
bedGraph2Cluster("bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph", "bedgraph/INPUT.WT.bedgraph,bedgraph/INPUT.dCDK.bedgraph", "bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph", "bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph,bedgraph/E2F1.bedgraph,bedgraph/CTCF.shSCR.bedgraph,bedgraph/c-Jun.shSCR.bedgraph", "output", "bed/hg19.200bp.bed", "2", "true", "8", "cosine", "1", "../")
```

## without license for MATLAB
[![PyPI version](https://badge.fury.io/py/run_matlab.svg)](https://badge.fury.io/py/run_matlab)
```shell
pip3 install run_matlab
run_matlab install -v R2018b -r 9.5
git clone https://github.com/hanjunlee21/bedGraph2Cluster
run_matlab run -v R2018b -r 9.5 bedGraph2Cluster bedGraph2Cluster bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph bedgraph/INPUT.WT.bedgraph,bedgraph/INPUT.dCDK.bedgraph bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph,bedgraph/E2F1.bedgraph,bedgraph/CTCF.shSCR.bedgraph,bedgraph/c-Jun.shSCR.bedgraph output bed/hg19.200bp.bed 2 true 8 cosine 1 $PWD
```

## BED_Bin

**hg19**: https://drive.google.com/uc?id=1chMyuUAK3rycgAbj1P-FnWhJkW7K0q0H&export=download&confirm=t

**hg38**: https://drive.google.com/uc?id=1kfVOx78xTBoJHC7EL_HnkCFiWxdrNDal&export=download&confirm=t


<br/>


# Examplary output
## All output files
![outputs](https://user-images.githubusercontent.com/67846757/187093037-7e33c73b-ceea-483e-a71c-79e791321718.png)

## Clustering heatmap
*k*=3             |  *k*=8
:-------------------------:|:-------------------------:
![fig1](https://user-images.githubusercontent.com/67846757/187092498-a6b716a2-1f01-4fc6-b437-298b60833b8f.png) | ![fig2](https://user-images.githubusercontent.com/67846757/187092500-4181af41-6641-4901-b9a9-eb26ec8345ad.png)


