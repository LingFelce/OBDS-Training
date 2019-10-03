"""
15_ChIP-seq workflow - 3rd October 2019

ChIP-seq analysis: identify TF binding sites (promoters/enhancers, target genes, sequence motifs, interacting partners, 
differential binding), chromatin state (active, poised, repressed)

Sources of bias:
- starting material - need many cells for robust results, consistent between replicates
- fragmentation of library (sonication - open chromatin shears more easily; enzymatic - sequence bias)
- antibody specificity. Use cocktail of antibodies - target different epitopes, avoid epitope masking. Use polyclonal Ab
- batch effects. 
ChIP-seq controls:
- mock (just Ab) - very little material - worth sequencing?
- no Ab (chromatin but IP with IgG) - very little material - worth sequencing?
- input control (chromatin not IP) - fragmentation bias
- spike in external control eg drosophila DNA
Sequencing considerations:
- single/paired end - paired end better for removal of duplicates
- read length long enough for unique mapping
- sequencing depth - greater for input (2 x IP), greater for broad peaks, Encode guidelines.
Mapping QC:
- > 90% mapping
- read filtering - mitochondrial reads, duplicates, multi-mapping reads, reads not properly paired
Encode metrics - PCR bottlenecking coefficients, non-redundant fraction
Peak calling:
- pool sample - better coverage
- caveats - samples should be of similar depth, otherwise one sample might influence peak calling more than others (down-sample?)
Peak calling tools - MACS2 (most widely used), SICER, LanceOtron, Homer, SPP, PeakRanger
Viewing peaks (in IGV) - different file formats:
- Bam (check raw reads, memory intensive) 
- Bedcoverage (breaks genome into regions and gives coverage score, make with Bedtools, plain text, less memory intensive)
- Bigwig (like bedcoverage, but compressed and indexed; faster, less memory intensive)
Peak calling QC:
- reproducibility (overlap replicates Bedtools, keep peaks 2+ samples)
- irreproducible discovery rate (Encode project), rank peaks on p-value, identify reproducible peaks on correlation of ranks
- black-listed regions (sticky regions that always appear in ChIP-seq)
Peak annotation:
- genomic context (promoter, enhancer etc) - GREAT online tool
- target genes (nearest neighbour) - GREAT
Differential binding (2 different conditions)
- featureCounts (count reads under peaks in each condition)
- DESeq2/edgeR (treat like RNA-seq data, must have >= 3 replicates)
Motif analysis: MEME suit, Homer motif analysis

"""
