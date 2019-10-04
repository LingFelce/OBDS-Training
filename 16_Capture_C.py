"""
16_Capture-C 4th October 2019

Chromosome conformation capture techniques (often abbreviated to 3C technologies or 3C-based methods) are a set of molecular biology 
methods used to analyze the spatial organization of chromatin in a cell. 
These methods quantify the number of interactions between genomic loci that are nearby in 3-D space, but may be separated by many 
nucleotides in the linear genome. Such interactions may result from biological functions, such as promoter-enhancer interactions, 
or from random polymer looping, where undirected physical motion of chromatin causes loci to collide.
Interaction frequencies may be analyzed directly, or they may be converted to distances and used to reconstruct 3-D structures.

First, the cell genomes are cross-linked with formaldehyde, which introduces bonds that "freeze" interactions between genomic loci. 
Treatment of cells with 1-3% formaldehyde, for 10-30min at room temperature is most common, however, standardization for preventing high 
protein-DNA cross linking is necessary, as this may negatively affect the efficiency of restriction digestion in the subsequent step. 
The genome is then cut into fragments with a restriction endonuclease. The size of restriction fragments determines the resolution of 
interaction mapping. Restriction enzymes (REs) that make cuts on 6bp recognition sequences, such as EcoR1 or HindIII, are used for this 
purpose, as they cut the genome once every 4000bp, giving ~ 1 million fragments in the human genome. For more precise interaction 
mapping, a 4bp recognizing RE may also be used. The next step is, random ligation. This takes place at low DNA concentrations in the 
presence of T4 DNA ligase, such that ligation between cross-linked interacting fragments is favored over ligation between fragments 
that are not cross-linked. 
Subsequently, interacting loci are quantified by amplifying ligated junctions by PCR methods

Sequence capture-based methods
A number of methods use oligonucleotide capture to enrich 3C and Hi-C libraries for specific loci of interest.
These methods include Capture-C, NG Capture-C, Capture-3C, and Capture Hi-C.
These methods are able to produce higher resolution and sensitivity than 4C based methods.

Capture-C similar to 4C:
Chromosome conformation capture-on-chip (4C) captures interactions between one locus and all other genomic loci. 
It involves a second ligation step, to create self-circularized DNA fragments, which are used to perform inverse PCR. 
Inverse PCR allows the known sequence to be used to amplify the unknown sequence ligated to it.
In contrast to 3C and 5C, the 4C technique does not require the prior knowledge of both interacting chromosomal regions. 
Results obtained using 4C are highly reproducible with most of the interactions that are detected between regions proximal to one another. 
On a single microarray, approximately a million interactions can be analyzed.

Use biotinylated probes to enrich for captured sequences, don't sequence whole library, just sites of interest.

See Jelena's notes:
http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/ppMan/CCseqBasic/2_workflow/index.html
http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/ppMan/CCseqBasic/DOCS/readsFragments_multipage.pdf
http://userweb.molbiol.ox.ac.uk/public/telenius/captureManual/UserManualforCaptureCanalysis.pdf
"""
#Capture-C pipeline exercise - files required: read1.fastq and read2.fastq, fragment.txt Hba-1 mouse globin locus

#fastqc code as before

#trimming using trim-galore - need to use collate decorator as need to put in reads as pair for trimming
#flash for fast length adjustment for short reads

