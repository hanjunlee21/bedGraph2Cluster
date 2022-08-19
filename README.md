# bedGraph2Cluster
MATLAB Function for the k-Means Clustering of ChIP-seq bedGraph Data

## Required arguments
    bedGraphs_Target: path(s) to target bedGraph files (included during clustering)
        possible values: string,
                         character array,
                         matrix of strings,
                         matrix of character arrays
    bedGraphs_Nontarget: path(s) to nontarget bedGraph files (excluded during clustering)
        possible values: string,
                         character array,
                         matrix of strings,
                         matrix of character arrays
    bedGraphs_Control: path(s) to control bedGraph files
        possible values: string,
                         character array,
                         matrix of strings,
                         matrix of character arrays
    BED_Bin: path to bin BED file
        possible values: string,
                         character array
    Outdir: path to output directory
        possible values: string,
                         character array
    k: number of clusters
        possible values: integer

## Optional arguments
    QNorm: whether to perform the QNorm normalization method
        possible values: "true",
                         "false",
                         'true',
                         'false',
                         true,
                         false
        caution: leaving QNorm empty would run bedGraph2Cluster
                 with QNorm normalization by default
        caution: opting for false would run bedGraph2Cluster
                 with CPM normalization method instead
