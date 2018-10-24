- General steps in RNA-Seq processing
- Designed for longRNA-Seq (>200bp)
- A lot of steps and settings are used as in https://github.com/ENCODE-DCC/long-rna-seq-pipeline ( https://www.encodeproject.org/pipelines/ ) with few modification
- Other possible source (not used) is https://github.com/ewels/NGI-RNAseq/blob/master/main.nf
- Quality steps/checks are used as in from https://github.com/SciLifeLab/NGI-RNAseq with few modification
- Should handle SE and PE as well 
- alignment NOT designed for QuantSeq and other "specialized" protocols" (but most of the QC should work)
- preprocessing, alignment, gene counts (probably) and QC NOT designed for shortRNA-Seq

TODO
- add general samples heatmaps
- add R scripts for DE
- add eXpress analysis for low count isoform expression and for another estimation type (https://github.com/COMBINE-lab/salmon/issues/107, http://www.nature.com/nmeth/journal/v10/n1/fig_tab/nmeth.2251_F3.html)
- add automatic analysis of dupRadar genes covered solely by multimapped reads (make table with common gene names)
- add StringTie for genome based transcriptome assembly https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual. Might need HISAT2 for alignment?
        script:
        def st_direction = ''
        if (forward_stranded && !unstranded){
            st_direction = "--fr"
        } else if (reverse_stranded && !unstranded){
            st_direction = "--rf"
        }
        """
        stringtie $bam_stringtieFPKM \\
            $st_direction \\
            -o ${bam_stringtieFPKM.baseName}_transcripts.gtf \\
            -v \\
            -G $gtf \\
            -A ${bam_stringtieFPKM.baseName}.gene_abund.txt \\
            -C ${bam_stringtieFPKM}.cov_refs.gtf \\
            -e \\
            -b ${bam_stringtieFPKM.baseName}_ballgown
- check Salmon and/or Kallisto quantification
