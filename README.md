# PROMPTs-Finder
The pipeline of finding PROMoter uPstream Transcripts(PROMPTs)
Identification of PROMPTs
First, we set background area as intergenic regions which are 10kb away from annotated genes. 
To generate the empirical distribution of ChrRNA-seq background signals, we randomly selected 230,000 boxes (200bp each) within background regions and calculated the ChrRNA-seq density of each box, resulting in an empirical distribution function. Next, we used sliding window (200bp in length, 10bp steps) to scan across the genome.
Each sliding window have values, fold changes over average of background areas and the probability to be different from background based on our empirical distribution function (adjusted by false discovery rate (FDR)).
Windows with FDR > 0.05 were removed. Remaining Windows in upstream antisense region of active promoters were merged if gap between windows less than 400bp. 
The merged areas were identified as PROMPTs.
PROMTPs added into hg38 ncbiRefSeq gtf file to provide features, and features expression quantification was performed with featureCounts tool of Subread (v2.0.1).
Differential expression analysis was performed using DESeq2 R package (v1.34.0). 
PROMPTs further selected with FDR > 0.05, FC > 2 (INTS11KO over WT) and start within 500bp upstream transcript start site.
