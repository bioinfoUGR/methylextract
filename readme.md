MethylExtract
-------------

Whole genome methylation profiling at a single cytosine resolution is now feasible due to the advent of high-throughput sequencing techniques together with bisulfite treatment of the DNA. To obtain the methylation value of each individual cytosine, the bisulfite treated sequence reads are first aligned to a reference genome, profiling afterwards the methylation levels from the alignments. A huge effort has been made to fast and correctly align the reads and many different algorithms and programs do exist. However, the second step is likewise crucial and non-trivial, but much less attention was paid to the final inference of the methylation states. Important error sources do exist like sequencing errors, bisulfite failure, clonal reads and single nucleotide variants.

We developed MethylExtract, a user friendly tool to generate i) high quality, whole genome methylation maps and ii) to detect sequence variation within the same sample preparation. The program is implemented into a single script and takes into account all major error sources. MethylExtract detects variation (SNVs – Single Nucleotide Variation) in a similar way than VarScan, a very sensitive method extensively used in SNP and genotype calling based on non-bisulfite treated reads. The usefulness of MethylExtract is shown by means of extensive benchmarking based on artificial BS reads and a comparison to a recently published method, Bis-SNP.

MethylExtract is able to detect SNVs within High-Throughput Sequencing experiments of bisulfite treated DNA at the same time as it generates high quality methylation maps. This simultaneous detection of DNA methylation and sequence variation is crucial for many downstream analyses, for example when deciphering the impact of Single Nucleotide Variants on differential methylation. An exclusive feature of MethylExtract, in comparison with existing software, is the possibility to assess in a statistical way the bisulfite failure. 

**Examples:** http://bioinfo2.ugr.es/MethylExtract/downloads/Example.tgz

**Artificial BS datasets:** http://bioinfo2.ugr.es/MethylExtract/downloads/ArtificialDatasets.tgz

**Reference**

Barturen G, Rueda A, Oliver JL and Hackenberg M (2014) MethylExtract: High-Quality methylation maps and SNV calling from whole genome bisulfite sequencing data [v2; ref status: indexed, http://f1000r.es/301] F1000Research 2014, 2:217 (doi:10.12688/f1000research.2-217.v2)



    
