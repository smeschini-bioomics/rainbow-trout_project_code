conda activate qiime2-amplicon-2024.10

qiime metadata tabulate \
 --m-input-file metadata.tsv \
 --o-visualization metadata.qzv

#Data importing in qiime environnement
qiime tools import \                                            
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path STPN2309_16S_samples_used_for_analysis \
 --input-format CasavaOneEightLanelessPerSampleDirFmt \
 --output-path STPN2309_demux-paired-end.qza
qiime demux summarize \
 --i-data STPN2309_demux-paired-end.qza \
 --o-visualization STPN2309_demux-paired-end.qzv

#trimming primers from sequences
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences STPN2309_demux-paired-end.qza \
 --p-front-f CCTACGGGNGGCWGCAG \
 --p-front-r GACTACHVGGGTATCTAATCC \
 --p-overlap 15 \
 --p-cores 12 \
 --p-match-adapter-wildcards \
 --o-trimmed-sequences STPN2309_primers_trimmed-seqs.qza \
 --verbose
qiime demux summarize \
 --i-data STPN2309_primers_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_trimmed-seqs.qza \
  --output-path exported_fastqs_primers-trimmed

#trimming adapters from TruqSeq kits 
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences STPN2309_primers_trimmed-seqs.qza \
 --p-adapter-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 --p-adapter-r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 --p-overlap 25 \
 --p-match-adapter-wildcards \
 --p-cores 12 \
 --o-trimmed-sequences STPN2309_primers_adapter_trimmed-seqs.qza \
 --verbose
qiime demux summarize \
 --i-data STPN2309_primers_adapter_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_adapter_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_adapter_trimmed-seqs.qza \
  --output-path exported_fastqs_primers_adapters-trimmed

#trimming poly-G and poly-A tails from sequences
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences STPN2309_primers_adapter_trimmed-seqs.qza \
  --p-adapter-f GGGGGGGGGGGGGGGGGG \
  --p-adapter-r GGGGGGGGGGGGGGGGGG \
  --p-adapter-f AAAAAAAAAAAAAAAAAA \
  --p-adapter-r AAAAAAAAAAAAAAAAAA \
  --p-overlap 10 \
  --p-match-adapter-wildcards \
  --p-cores 12 \
  --o-trimmed-sequences STPN2309_primers_adapters_poly_tail_trimmed-seqs.qza \
  --verbose
qiime demux summarize \
 --i-data STPN2309_primers_adapters_poly_tail_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_adapters_poly_tail_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_adapters_poly_tail_trimmed-seqs.qza \
  --output-path exported_fastqs_primers_adapters_poly_tail-trimmed

#keep the sequences that have a minimum lenght of 180 bp
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences STPN2309_primers_adapters_poly_tail_trimmed-seqs.qza \
 --p-minimum-length 180 \
 --p-cores 12 \
 --o-trimmed-sequences STPN2309_primers_adapters_poly_tail_180ml_trimmed-seqs.qza \
 --verbose
qiime demux summarize \
 --i-data STPN2309_primers_adapters_poly_tail_180ml_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_adapters_poly_tail_180ml_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_adapters_poly_tail_180ml_trimmed-seqs.qza \
  --output-path exported_fastqs_primers_adapters_poly_tail_180_ml-trimmed

#generating a FastaQC (for each sample) and multiQC (global R1 and R2) quality report to evaluate trimming 
cd PY
python python_script_fastaQC_report.py
python python_script_multiQC_report.py

#after multiQC report, it stills some adapter sequences (<8%) in some reads that have not been removed because of minimum of 25bp of overlapping in previous command line. So, a less stringent trimming will be applied (10 of overlapping)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences STPN2309_primers_adapters_poly_tail_180ml_trimmed-seqs.qza \
 --p-adapter-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
 --p-adapter-r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
 --p-overlap 10 \
 --p-match-adapter-wildcards \
 --p-cores 12 \
 --p-minimum-length 180 \
 --o-trimmed-sequences STPN2309_primers_adapters_poly_tail_180ml_multiqc_trimmed-seqs.qza \
 --verbose
qiime demux summarize \
 --i-data STPN2309_primers_adapters_poly_tail_180ml_multiqc_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_adapters_poly_tail_180ml_multiqc_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_adapters_poly_tail_180ml_multiqc_trimmed-seqs.qza \
  --output-path exported_fastqs_primers_adapters_poly_tail_180_ml_multiqc-trimmed
  
#generating again a FastaQC (for each sample) and multiQC (global R1 and R2) quality report to evaluate trimming (remember to change the dir input name for the new reports!)
cd PY
python python_script_fastaQC_report.py
python python_script_multiQC_report.py

# trim low-quality bases from the 3' end based on a Phred score (below 20)
qiime cutadapt trim-paired \
 --i-demultiplexed-sequences STPN2309_primers_adapters_poly_tail_180ml_multiqc_trimmed-seqs.qza \
 --p-quality-cutoff-3end 20 \
 --o-trimmed-sequences STPN2309_primers_adapters_poly_tail_180ml_multiqc_20quality_trimmed-seqs.qza \
 --p-cores 12 \
 --verbose
qiime demux summarize \
 --i-data STPN2309_primers_adapters_poly_tail_180ml_multiqc_20quality_trimmed-seqs.qza \
 --o-visualization STPN2309_primers_adapters_poly_tail_180ml_multiqc_20quality_trimmed-seqs.qzv
qiime tools export \
  --input-path STPN2309_primers_adapters_poly_tail_180ml_multiqc_20quality_trimmed-seqs.qza \
  --output-path exported_fastqs_primers_adapters_poly_tail_180_ml_20quality-trimmed

#denoise-paired: Denoise and dereplicate paired-end sequences (don't select multithreading or an eccessive number of cores because crash memory issue with related errors can happens! "'names' attribute [96] must be the same length as the vector [92]")
  --i-demultiplexed-seqs STPN2309_primers_adapters_poly_tail_180ml_multiqc_20quality_trimmed-seqs.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-table STPN2309_table_after_dada2.qza \
  --o-representative-sequences STPN2309_rep-seqs_after_dada2.qza \
  --o-denoising-stats STPN2309_denoising-stats.qza \
  --verbose
qiime feature-table summarize \
  --i-table STPN2309_table_after_dada2.qza \
  --o-visualization STPN2309_table_after_dada2.qzv \
  --m-sample-metadata-file STPN2309_16S_metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data STPN2309_rep-seqs_after_dada2.qza \
  --o-visualization STPN2309_rep-seqs_after_dada2.qzv
qiime metadata tabulate \
  --m-input-file STPN2309_denoising-stats.qza \
  --o-visualization STPN2309_denoising-stats.qzv

### Processing and filtering SILVA database with RESCRIPt using 'extract-seq-segment' (-> pool expansion steps)
qiime rescript get-silva-data \
 --p-version '138.2' \
 --p-target 'SSURef_NR99' \
 --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
 --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza
qiime rescript reverse-transcribe \
 --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
 --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza
qiime rescript cull-seqs \
 --i-sequences silva-138.2-ssu-nr99-seqs.qza \
 --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza
qiime rescript filter-seqs-length-by-taxon \
 --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
 --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
 --p-labels Archaea Bacteria Eukaryota \
 --p-min-lens 900 1200 1400 \
 --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
 --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza
qiime rescript dereplicate \
 --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
 --i-taxa silva-138.2-ssu-nr99-tax.qza \
 --p-mode 'uniq' \
 --p-threads 8 \
 --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
 --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza
# Generate an initial reference pool of sequence segments using primer-pair search
qiime feature-classifier extract-reads \
 --i-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
 --p-f-primer CCTACGGGNGGCWGCAG \
 --p-r-primer GACTACHVGGGTATCTAATCC \
 --p-n-jobs 8 \
 --p-read-orientation 'forward' \
 --o-reads silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq.qza
qiime rescript dereplicate \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq.qza \
 --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
 --p-mode 'uniq' \
 --p-threads 8 \
 --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep.qza \
 --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-segments-uniq.qza
qiime rescript cull-seqs \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep.qza  \
 --p-n-jobs 8 \
 --p-num-degenerates 1 \
 --p-homopolymer-length 8 \
 --o-clean-sequences silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep-cull.qza
qiime feature-table tabulate-seqs \
 --i-data silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep-cull.qza \
 --o-visualization silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep-cull.qzv
# First iteration of query sequence segment extraction (-> pool expansion step 1)
qiime rescript extract-seq-segments \
 --i-input-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
 --i-reference-segment-sequences silva-138.2-ssu-nr99-seqs-V3-V4-mode-uniq-segment-derep-cull.qza \
 --p-perc-identity 0.7 \
 --p-min-seq-len 10 \
 --p-threads 8 \
 --o-extracted-sequence-segments silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-01.qza \
 --o-unmatched-sequences silva-138.2-ssu-nr99-seqs-V3-V4-unmatched-sequences-01.qza \
 --verbose
qiime rescript dereplicate \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-01.qza \
 --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
 --p-mode 'uniq' \
 --p-threads 8 \
 --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-01.qza \
 --o-dereplicated-taxa silva-138.2-ssu-nr99-extracted-tax-segments-derep-01.qza
qiime rescript cull-seqs \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-01.qza \
 --p-n-jobs 8 \
 --p-num-degenerates 1 \
 --p-homopolymer-length 8 \
 --o-clean-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-01.qza
qiime feature-table tabulate-seqs \
 --i-data silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-01.qza \
 --o-visualization silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-01.qzv
# The second iteration (-> pool expansion step 2)
qiime rescript extract-seq-segments \
 --i-input-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
 --i-reference-segment-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-01.qza \
 --p-perc-identity 0.7 \
 --p-min-seq-len 10 \
 --p-threads 8 \
 --o-extracted-sequence-segments silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-02.qza \
 --o-unmatched-sequences silva-138.2-ssu-nr99-seqs-V3-V4-unmatched-sequences-02.qza \
 --verbose
qiime rescript dereplicate \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-02.qza \
 --i-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza \
 --p-mode 'uniq' \
 --p-threads 8 \
 --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-02.qza \
 --o-dereplicated-taxa silva-138.2-ssu-nr99-extracted-tax-segments-derep-02.qza
qiime rescript cull-seqs \
 --i-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-02.qza \
 --p-n-jobs 8 \
 --p-num-degenerates 1 \
 --p-homopolymer-length 8 \
 --o-clean-sequences silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-02.qza
qiime feature-table tabulate-seqs \
 --i-data silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-02.qza \
 --o-visualization silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-02.qzv
qiime rescript filter-taxa \
 --i-taxonomy silva-138.2-ssu-nr99-extracted-tax-segments-derep-02.qza \
 --m-ids-to-keep-file silva-138.2-ssu-nr99-seqs-V3-V4-extracted-seq-segments-derep-cull-02.qza \
 --o-filtered-taxonomy silva-138.2-ssu-nr99-extracted-tax-segments-derep-cull-keep-02.qza
#Build and evaluate our classifier
qiime rescript evaluate-fit-classifier \
 --i-sequences 2024.09-extracted-seq-segments-derep-cull-02.qza  \
 --i-taxonomy 2024.09-extracted-tax-segments-derep-cull-keep-02.qza \
 --p-n-jobs 2 \
 --o-classifier gg2024.09-V3-V4-classifier.qza \
 --o-evaluation gg2024.09-V3-V4-classifier-evaluation.qzv \
 --o-observed-taxonomy gg2024.09-extracted-tax-refseqs-predicted-taxonomy.qza
qiime rescript evaluate-taxonomy \
 --i-taxonomies 2024.09-extracted-tax-segments-derep-cull-keep-02.qza gg2024.09-extracted-tax-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats gg2024.09-V3-V4-taxonomy-evaluation.qzv

# identify decontam in digesta samples
qiime feature-table filter-samples \
  --i-table STPN2309_table_after_dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[sampletype] = 'Digesta'" \
  --o-filtered-table STPN2309_table_after_dada2-digesta.qza
qiime quality-control decontam-identify  \
--i-table STPN2309_table_after_dada2-digesta.qza  \
--m-metadata-file metadata.tsv \
--p-method frequency  \
--p-freq-concentration-column Concentration  \
--o-decontam-scores digesta-freq_decontam_scores.qza
# inspect output
qiime quality-control decontam-score-viz  \
--i-decontam-scores digesta-freq_decontam_scores.qza  \
--i-table STPN2309_table_after_dada2-digesta.qza --i-rep-seqs STPN2309_rep-seqs_after_dada2.qza  \
--p-threshold 0.1  \
--p-no-weighted  \
--p-bin-size 0.05  \
--o-visualization digesta-decontam_score.qzv
#remove contaminants from digesta feature table
qiime feature-table filter-features  \
--i-table STPN2309_table_after_dada2-digesta.qza  \
--m-metadata-file digesta-freq_decontam_scores.qza  \
--p-where '[p]>0.1 OR [p] IS NULL'  \
--o-filtered-table STPN2309_decontam-filtered-table-digesta.qza
#remove feature from rep-seq table
qiime feature-table filter-seqs  \
--i-data STPN2309_rep-seqs_after_dada2.qza  \
--i-table STPN2309_decontam-filtered-table-digesta.qza  \
--o-filtered-data STPN2309_decontam-filtered_rep-seq-digesta.qza

# identify decontam in mucus samples
qiime feature-table filter-samples \
  --i-table STPN2309_table_after_dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[sampletype] = 'Mucus'" \
  --o-filtered-table STPN2309_table_after_dada2-mucus.qza
qiime quality-control decontam-identify  \
--i-table STPN2309_table_after_dada2-mucus.qza  \
--m-metadata-file metadata.tsv \
--p-method frequency  \
--p-freq-concentration-column Concentration  \
--o-decontam-scores mucus-freq_decontam_scores.qza
# inspect output
qiime quality-control decontam-score-viz  \
--i-decontam-scores mucus-freq_decontam_scores.qza  \
--i-table STPN2309_table_after_dada2-mucus.qza --i-rep-seqs STPN2309_rep-seqs_after_dada2.qza  \
--p-threshold 0.1  \
--p-no-weighted  \
--p-bin-size 0.05  \
--o-visualization mucus-decontam_score.qzv
#remove contaminants from digesta feature table
qiime feature-table filter-features  \
--i-table STPN2309_table_after_dada2-mucus.qza  \
--m-metadata-file mucus-freq_decontam_scores.qza  \
--p-where '[p]>0.1 OR [p] IS NULL'  \
--o-filtered-table STPN2309_decontam-filtered-table-mucus.qza
#remove feature from rep-seq table
qiime feature-table filter-seqs  \
--i-data STPN2309_rep-seqs_after_dada2.qza  \
--i-table STPN2309_decontam-filtered-table-mucus.qza  \
--o-filtered-data STPN2309_decontam-filtered_rep-seq-mucus.qza

# identify decontam in feed samples
qiime feature-table filter-samples \
  --i-table STPN2309_table_after_dada2.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[sampletype] = 'Feed'" \
  --o-filtered-table STPN2309_table_after_dada2-feed.qza

qiime quality-control decontam-identify  \
--i-table STPN2309_table_after_dada2-feed.qza  \
--m-metadata-file metadata.tsv \
--p-method frequency  \
--p-freq-concentration-column Concentration  \
--o-decontam-scores feed-freq_decontam_scores.qza
# inspect output
qiime quality-control decontam-score-viz  \
--i-decontam-scores feed-freq_decontam_scores.qza  \
--i-table STPN2309_table_after_dada2-feed.qza --i-rep-seqs STPN2309_rep-seqs_after_dada2.qza  \
--p-threshold 0.1  \
--p-no-weighted  \
--p-bin-size 0.05  \
--o-visualization feed-decontam_score.qzv
#remove contaminants from digesta feature table
qiime feature-table filter-features  \
--i-table STPN2309_table_after_dada2-feed.qza  \
--m-metadata-file feed-freq_decontam_scores.qza  \
--p-where '[p]>0.1 OR [p] IS NULL'  \
--o-filtered-table STPN2309_decontam-filtered-table-feed.qza
#remove feature from rep-seq table
qiime feature-table filter-seqs  \
--i-data STPN2309_rep-seqs_after_dada2.qza  \
--i-table STPN2309_decontam-filtered-table-feed.qza  \
--o-filtered-data STPN2309_decontam-filtered_rep-seq-feed.qza

# merge tables
qiime feature-table merge \
--i-tables STPN2309_decontam-filtered-table-digesta.qza \
--i-tables STPN2309_decontam-filtered-table-mucus.qza \
--i-tables STPN2309_decontam-filtered-table-feed.qza \
--o-merged-table STPN2309_decontam-filtered-merged-table.qza
# merge rep seq
qiime feature-table merge-seqs \
--i-data STPN2309_decontam-filtered_rep-seq-digesta.qza \
--i-data STPN2309_decontam-filtered_rep-seq-mucus.qza \
--i-data STPN2309_decontam-filtered_rep-seq-feed.qza \
--o-merged-data STPN2309_decontam-filtered-merged-rep-seqs.qza




#Get the taxonomy
qiime feature-classifier classify-sklearn \
 --i-classifier silva-138.2-ssu-nr99-seqs-V3-V4-classifier.qza \
 --i-reads STPN2309_decontam-filtered-merged-rep-seqs.qza \
 --o-classification STPN2309-taxonomy-silva-138.2.qza
qiime metadata tabulate \
 --m-input-file STPN2309-taxonomy-silva-138.2.qza \
 --o-visualization STPN2309-taxonomy-silva-138.2.qzv

#filter table
qiime feature-table filter-features \
  --i-table STPN2309_decontam-filtered-merged-table.qza \
  --p-min-samples 2 \
  --o-filtered-table STPN2309_decontam-filtered-merged-table-freq1.qza

qiime taxa filter-table \
  --i-table STPN2309_decontam-filtered-merged-table-freq1.qza \
  --i-taxonomy STPN2309-taxonomy-silva-138.2.qza \
  --p-include p__ \
  --o-filtered-table STPN2309_decontam-filtered-merged-table-freq1-with-phyla.qza

qiime feature-table summarize \
 --i-table STPN2309_decontam-filtered-merged-table-freq1-with-phyla.qza \
 --o-visualization STPN2309_decontam-filtered-merged-table-freq1-with-phyla.qzv

# synchronized rep seq
qiime feature-table filter-seqs  \
 --i-data STPN2309_decontam-filtered-merged-rep-seqs.qza  \
 --i-table STPN2309_decontam-filtered-merged-table-freq1-with-phyla.qza  \
 --o-filtered-data STPN2309_decontam-filtered-merged-rep-seqs-freq1-with-phyla.qza

#Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
 --i-sequences STPN2309_decontam-filtered-merged-rep-seqs-freq1-with-phyla.qza \
 --o-alignment STPN2309-silva_nb-aligned.qza \
 --o-masked-alignment STPN2309-silva_nb-masked.qza \
 --o-tree STPN2309-silva_nb-unrooted-tree.qza \
 --o-rooted-tree STPN2309-silva_nb-rooted-tree.qza

