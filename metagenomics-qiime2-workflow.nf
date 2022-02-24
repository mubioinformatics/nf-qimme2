#!/usr/bin/env nextflow

params.target_dir = ""
params.marker = "16S_V4_bacteria"
//params.classifier_qza = "V416S_classifier.current"
params.cpus = 30
params.poolname = "Test"
params.references = "References.txt"

//Give flexibility to decide on classifier based on command-line argument
if(params.marker == "16S_V4_bacteria"){
    //CLASSIFIER = '/data/references/classifiers/16S_V4/16S_V4_classifier.qza'
    //CLASSIFIER_INFO = "Green Genes, v13_8, 99% for library preps of PE250"
    CLASSIFIER = '/data/references/classifiers/16S_V4/16S_V4_silva_PE250_classifier.qza'
    CLASSIFIER_INFO = "Silva v138, 99%, for library preps of PE250"
    FRONT_F = "GTGCCAGCMGCCGCGGTAA"
    ADAPTER_F = "ATTAGAWACCCBDGTAGTCC"
    FRONT_R = "GGACTACHVGGGTWTCTAAT"
    ADAPTER_R = "TTACCGCGGCKGCTGGCAC"
}
else{
    CLASSIFIER = '/data/references/classifiers/its/ITS_Unite99_classifier.qza'
    CLASSIFIER_INFO = "Unite 10.05.2021, 99% (no primer removal)"
    FRONT_F = "GTCCCTGCCCTTTGTACACA" 
    ADAPTER_F = "CGATGAAGAACGCAGCGAAA" //rev comp of front-r
    FRONT_R = "TTTCGCTGCGTTCTTCATCG"
    ADAPTER_R = "TGTGTACAAAGGGCAGGGAC" //rev comp of front-f
    //ADAPTER_R = "GATCCYTCCGCAGGTAG" //old
}

params.manifest = ""
params.metadata = ""
qiime="PATH=/cluster/software/qiime2-2021.8/bin:$PATH qiime"
biom="PATH=/cluster/software/qiime2-2021.8/bin:$PATH biom"

process demux {
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else null
            }

    input:
    path manifest from params.manifest

    output:
    path "paired-end-demux.qza" into demux_art
    path "paired-end-demux.qzv" into demux_viz

    script:
    """
    export MPLCONFIGDIR=${workDir}
    $qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path $manifest \
        --output-path paired-end-demux.qza \
        --input-format PairedEndFastqManifestPhred33

    $qiime demux summarize \
        --i-data paired-end-demux.qza \
        --o-visualization paired-end-demux.qzv
    """
}

process cutadapt {
    //publishDir "${baseDir}/${params.poolname}/Demux", mode: 'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else if (fn.indexOf("demux-s-primers.qza") > 0) null
            else "Demux/$fn"
        }

    input:
    path input_qza from demux_art

    output:
    path "paired-end-demux-s-primers.qza" into cutadapt_art
    path "paired-end-demux-s-primers.qzv" into cutadapt_viz

    script:
    """
    
    $qiime cutadapt trim-paired \
        --i-demultiplexed-sequences  $input_qza \
        --p-cores ${params.cpus} \
        --p-front-f ${FRONT_F} \
        --p-adapter-f ${ADAPTER_F} \
        --p-front-r ${FRONT_R} \
        --p-adapter-r ${ADAPTER_R} \
        --p-no-indels \
        --p-times 1 \
        --o-trimmed-sequences paired-end-demux-s-primers.qza \
        --p-discard-untrimmed \
        --verbose

    # Cutadapt result visualization
    $qiime demux summarize \
        --i-data paired-end-demux-s-primers.qza \
        --o-visualization paired-end-demux-s-primers.qzv
    """
}

process dada2 {
    //publishDir "${baseDir}/${params.poolname}/Dada2", mode: 'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else "Dada2/$fn"
        }

    input:
    path input_qza from cutadapt_art

    output:
    path "representative_sequences.qza" into dada2_repseqs_art, dada2_repseqs_art_biom, dada2_repseq_art_phylo
    path "representative_sequences.qzv" into dada2_repseq_viz
    path "dada2_feature_table.qza" optional true into dada2_feature_table_art
    path "table.qza" into dada2_table_art, \
                        dada2_table_art_forBiom, \
                        dada2_table_art_forConsolidate, \
                        dada2_table_art_forDiversity, \
                        dada2_table_for_citation
    path "Denoising_stats.qza" into denoising_stats_art
    path "table.qzv" into dada2_table_viz
    path "dna-sequences.fasta" optional true into dada2_dna_seq_ch
    path "accepted.seq.lengths" optional true into dada2_acc_seqlength_ch

    shell:
    if( params.marker == '16S_V4_bacteria' )
    """
    $qiime dada2 denoise-paired \
        --i-demultiplexed-seqs !{input_qza} \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-n-threads !{params.cpus} \
        --p-chimera-method consensus \
        --o-table dada2_feature_table.qza \
        --output-dir Dada2 \
        --o-representative-sequences representative_sequences.qza \
        --o-denoising-stats Denoising_stats.qza \
        --verbose

    mv representative_sequences.qza original_representative_sequences.qza
    #mv dada2_feature_table.qza original_table.qza

    # export sequence data
    $qiime tools export \
        --input-path original_representative_sequences.qza \
        --output-path ./dada2_tmp

    # create set of sequences with desired length range
    cat dada2_tmp/dna-sequences.fasta |\
    sed '/^>.*\$/N; s/\\n/\\t/' |\
    awk '{
            if(length(\$2)>=249 && length(\$2)<=257) print ">"\$1"\\n"\$2
    }' | sed 's/>//' >  dna-sequences.fasta

    rm -fr dada2_tmp

    cat dna-sequences.fasta | grep '>' | sed 's/>//' |\
        awk 'BEGIN{print "feature-id"}{print \$1}' > accepted.seq.lengths

    $qiime tools import --type FeatureData[Sequence] \
        --input-path dna-sequences.fasta \
        --output-path representative_sequences.qza


    # With the list of seqids with correct seq lengths, trim table
    $qiime feature-table filter-features \
        --i-table dada2_feature_table.qza \
        --m-metadata-file accepted.seq.lengths \
        --o-filtered-table table.qza

    $qiime feature-table summarize \
        --i-table table.qza \
        --m-sample-metadata-file !{params.metadata} \
        --o-visualization table.qzv
        
    $qiime feature-table tabulate-seqs \
        --i-data representative_sequences.qza \
        --o-visualization representative_sequences.qzv
    """
    //removed: --output-dir Dada2 \ because nextflow handles this
    else
    """
    $qiime dada2 denoise-paired \
        --i-demultiplexed-seqs !{input_qza} \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-max-ee-f 2 \
        --p-max-ee-r 2 \
        --p-n-threads !{params.cpus} \
        --p-chimera-method consensus \
        --o-table table.qza \
        --output-dir Dada2 \
        --o-representative-sequences representative_sequences.qza \
        --o-denoising-stats Denoising_stats.qza \
        --verbose

    $qiime feature-table summarize \
        --i-table table.qza \
        --m-sample-metadata-file !{params.metadata} \
        --o-visualization table.qzv

    $qiime feature-table tabulate-seqs \
        --i-data representative_sequences.qza \
        --o-visualization representative_sequences.qzv
    """
}

process classify_and_plot {
    //publishDir "${baseDir}/${params.poolname}", mode: 'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else "Classification/$fn"
        }
    
    input:
    path input_repseq_qza from dada2_repseqs_art
    path input_table_qza from dada2_table_art

    output:
    path "Classifier_output.qza" into classifier_output_art,classifier_cite
    path "taxa-bar-plots.qzv" into taxa_barplot_viz

    script:
    """
    $qiime feature-classifier classify-sklearn \
        --i-reads ${input_repseq_qza} \
        --i-classifier $CLASSIFIER \
        --p-n-jobs ${params.cpus} \
        --o-classification Classifier_output.qza \
        --verbose

    $qiime taxa barplot \
        --i-table ${input_table_qza} \
        --i-taxonomy Classifier_output.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization taxa-bar-plots.qzv
    """
}

process biom_collate {
    publishDir "${baseDir}/${params.poolname}/", mode: 'copy'

    input:
    path table_qza from dada2_table_art_forBiom
    path classifier_output from classifier_output_art
    path repseq_qza from dada2_repseqs_art_biom
    path paired_end_demux_qzv from demux_viz
    path denoising_stats_art from denoising_stats_art

    output:
    path "Biom_dir" into biom_dir,biom_dir_rare,biom_dir_alpha,biom_make_csv

    script:
    """
    $qiime tools export --input-path ${table_qza} --output-path Biom_dir
    $qiime tools export --input-path ${classifier_output} --output-path Biom_dir

    $biom add-metadata \
        -i Biom_dir/feature-table.biom \
        -o Biom_dir/table_with_obs_metadata.biom \
        --observation-metadata-fp Biom_dir/taxonomy.tsv \
        --observation-header "Feature ID,Taxon" \
        --sc-separated taxonomy
    
    $biom convert -i Biom_dir/table_with_obs_metadata.biom \
        -o Biom_dir/annotated_table.csv \
        --to-tsv \
        --output-metadata-id=Taxon \
        --tsv-metadata-formatter=naive \
        --header-key=Taxon

    $biom summarize-table \
        -i Biom_dir/table_with_obs_metadata.biom \
        -o Biom_dir/Sample.summary.csv
    sed -i 's/\\:/\\t/' Biom_dir/Sample.summary.csv

    $biom summarize-table \
        -i Biom_dir/table_with_obs_metadata.biom \
        -o Biom_dir/ASV.summary.csv \
        --observations
    sed -i 's/\\:/\\t/' Biom_dir/ASV.summary.csv

    $qiime tools export --input-path ${repseq_qza} --output-path Biom_dir/rep-seqs
    
    cat Biom_dir/rep-seqs/dna-sequences.fasta |\
        awk '{if(\$1~/>/){
            if(NR>1)print ID,LEN"\\n"SEQ; 
            ID=\$1;
            LEN=0
            }
            else{
            LEN=length();
            SEQ=\$0
            }
        }END{print ID,LEN"\\n"SEQ}' >  Biom_dir/rep-seqs.csv
    
    rm -fr Biom_dir/rep-seqs

    $qiime tools export --input-path ${paired_end_demux_qzv} --output-path junk
    cat junk/per-sample-fastq-counts.tsv | \
        sed '1d;s/,/\\t/' | \
        sort -k1,1  > Biom_dir/raw-per-sample-fastq-counts.txt
    rm -fr junk

    cat Biom_dir/raw-per-sample-fastq-counts.txt |\
        awk 'BEGIN{print "sample-id\\tRawInput\\n#q2:types\\tnumeric"}{print \$0}' > Biom_dir/raw_cnts

    $qiime tools export --input-path ${denoising_stats_art} --output-path Denoising_stats
    mv Denoising_stats/stats.tsv Denoising_stats/stats.csv
    #dos2unix Denoising_stats/stats.csv

    cat  Biom_dir/Sample.summary.csv | \
    sed '1,/sample detail/d' | sort -k1,1 | sed 's/ //g' | \
        awk 'BEGIN{print "sample-id\\tLength_filtered:Counts\\n#q2:types\\tnumeric"}{
            print \$0
        }' | sed 's/,//g; s/\\.000\$//' > Biom_dir/Sample.summary.cnts

    # For Denoising_stats/stats.csv, Biom_dir/raw_cnts, and Biom_dir/Sample.summary.cnts
    # Capture the first two header lines, then sort the rest
    head -n2  Denoising_stats/stats.csv            > Denoising_stats/sorted_stats.csv
    tail -n+3 Denoising_stats/stats.csv    | sort >> Denoising_stats/sorted_stats.csv
    head -n2  Biom_dir/raw_cnts                    > Biom_dir/sorted_raw_cnts
    tail -n+3 Biom_dir/raw_cnts            | sort >> Biom_dir/sorted_raw_cnts
    head -n2  Biom_dir/Sample.summary.cnts         > Biom_dir/sorted_Sample.summary.cnts
    tail -n+3 Biom_dir/Sample.summary.cnts | sort >> Biom_dir/sorted_Sample.summary.cnts

    join -t \$'\\t' Biom_dir/sorted_raw_cnts  Denoising_stats/sorted_stats.csv | \
    sed '1s/input/From_Cutadapt/' | sed '1s/filtered/Dada2:filtered/' | \
    join -t \$'\\t' - Biom_dir/sorted_Sample.summary.cnts > Biom_dir/stats.csv
    #rm Biom_dir/Sample.summary.cnts Biom_dir/raw-per-sample-fastq-counts.txt Biom_dir/raw_cnts    
    """
}

process feature_grouping {
    //publishDir "${baseDir}/${params.poolname}", mode: 'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else "Grouped_by_Taxa/$fn"
        }

    input:
    path biom_dir from biom_dir
    path table_qza from dada2_table_art_forConsolidate

    output:
    path "grouped_by_taxa.qza" into grouped_by_taxa_art

    //May need a parameter to get all the things they'd like to group by
    script:
    """
    $qiime feature-table group \
        --i-table ${table_qza}  \
        --p-mode sum \
        --m-metadata-column Taxon \
        --m-metadata-file ${biom_dir}/taxonomy.tsv \
        --p-axis feature \
        --o-grouped-table grouped_by_taxa.qza
    """
}

process rarefaction {
    //publishDir "${baseDir}/${params.poolname}/rarefaction", mode: 'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else if(fn.indexOf(".errors") > 0) "$fn"
            else "Rarefaction/$fn"
        }

    input:
    path biom_dir from biom_dir_rare
    path grouped_taxa_qza from grouped_by_taxa_art

    output:
    path "rarefied-table.low.csv" into rare_low_depth_csv
    path "rarefied-table.high.csv" into rare_hi_depth_csv
    path "rarefaction.errors" optional true into rarefaction_errors

    //tail -n +3 Biom_dir/stats.csv| cut -f 12 | sort -rn | tail -n1
    script:
    """
    
    DEPTH1K=`cat ${biom_dir}/stats.csv | \
        awk 'BEGIN{THRESH=1000;RARE_DEPTH=1000000000}\
        NR>2{if(\$NF>THRESH && \$NF<RARE_DEPTH)RARE_DEPTH=\$NF}\
        END{print RARE_DEPTH-1}'`
    echo \$DEPTH1K

    $qiime feature-table rarefy \
        --i-table ${grouped_taxa_qza} \
        --p-sampling-depth \$DEPTH1K \
        --o-rarefied-table my-rarefied-table.1K.qza || echo "Failed to create low sampled rarefied table" >> rarefied-table.low.csv
    $qiime tools export --input-path my-rarefied-table.1K.qza --output-path . || echo "Creation of rarefied biom failed, probably failed to make low-sampled rarefied qza" >> rarefaction.errors
    $biom convert -i feature-table.biom -o rarefied-table.low.csv --to-tsv || echo "Creation of rarefied csv failed, probably failed to make low-sampled rarefied qza" >> rarefaction.errors

    DEPTH10K=`cat ${biom_dir}/stats.csv | \
	    awk 'BEGIN{THRESH=10000;RARE_DEPTH=1000000000}\
        NR>2{if(\$NF>THRESH && \$NF<RARE_DEPTH)RARE_DEPTH=\$NF}\
        END{print RARE_DEPTH-1}'`
    echo \$DEPTH10K

    $qiime feature-table rarefy \
            --i-table ${grouped_taxa_qza} \
            --p-sampling-depth \$DEPTH10K \
            --o-rarefied-table my-rarefied-table.10K.qza || echo "Failed to create high sampled rarefied table" >> rarefied-table.high.csv
    $qiime tools export --input-path my-rarefied-table.10K.qza --output-path . || echo "Creation of rarefied biom failed, probably failed to make high-sampled rarefied qza" >> rarefaction.errors
    $biom convert -i feature-table.biom -o rarefied-table.high.csv --to-tsv || echo "Creation of rarefied csv failed, probably failed to make high-sampled rarefied qza" >> rarefaction.errors
    """
}

process phylogeny {
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else "Tree/$fn"
        }

    input:
    path repseq_qza from dada2_repseq_art_phylo

    output:
    //path "aligned-rep-seqs.qza" into aligned_repseqs_art
    //path "masked-aligned-rep-seqs.qza" into masked_aligned_repseq_art
    //path "unrooted-tree.qza" into unrooted_tree_art
    path "rooted-tree.qza" into rooted_tree_art

    script:
    """
    $qiime alignment mafft \
        --p-n-threads ${params.cpus} \
        --i-sequences ${repseq_qza} \
        --o-alignment aligned-rep-seqs.qza

    #Filter the alignment to remove positions that are highly variable.
    $qiime alignment mask \
        --i-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza

    #Generate an un-rooted phylogenetic tree from the masked alignment
    $qiime phylogeny fasttree \
        --p-n-threads ${params.cpus} \
        --i-alignment masked-aligned-rep-seqs.qza \
        --o-tree unrooted-tree.qza

    #Place the root of the tree at the midpoint of the longest tip-to-tip distance in the unrooted tree
    $qiime phylogeny midpoint-root \
        --i-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
    """
}

process alpha_diversity {
    //publishDir "${baseDir}/${params.poolname}/", mode:'copy'
    publishDir "${baseDir}/${params.poolname}/",
        mode: 'copy',
            saveAs: {fn ->
            if (fn.indexOf(".qzv") > 0) "$fn"
            else if (fn.indexOf(".errors") >0) "$fn"
            else "$fn"
        }

    input:
    path biom_dir from biom_dir_alpha
    path rooted_tree from rooted_tree_art
    path table_qza from dada2_table_art_forDiversity

    output:
    path "alpha-rarefaction.qzv" into alpha_rarefaction_viz
    path "core-metrics-results" into core_metrics_outdir
    path "qiime_diversity.errors" into diversity_error_ch
    val 1 into completed_alpha_div

    script:
    """
    MAX_SMPL_CNT=`tail -n1 ${biom_dir}/Sample.summary.csv | cut -f2 | \
    sed 's/,//g' | cut -f1 -d '.'`    # Remove all commas

    $qiime diversity alpha-rarefaction \
        --i-table ${table_qza} \
        --i-phylogeny ${rooted_tree} \
        --p-max-depth \$MAX_SMPL_CNT \
        --m-metadata-file ${params.metadata} \
        --o-visualization alpha-rarefaction.qzv

    $qiime diversity core-metrics-phylogenetic \
        --i-table ${table_qza}  \
        --i-phylogeny ${rooted_tree}  \
        --p-n-jobs-or-threads ${params.cpus} \
        --p-sampling-depth 1000 \
        --m-metadata-file ${params.metadata} \
        --output-dir core-metrics-results

    $qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/faith-pd-group-significance.qzv 2>> qiime_diversity.errors ||\
            echo "Could not run alpha-group-significance on faith_pd_vector.qza, probably due to lack of categorical metadata column that meets visualizer's requirements." >> qiime_diversity.errors

    $qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/evenness_vector.qza \
        --m-metadata-file ${params.metadata}\
        --o-visualization core-metrics-results/evenness-group-significance.qzv 2>> qiime_diversity.errors ||\
            echo "Could not run alpha-group-significance on evenness_vector.qza, probably due to a lack of categorical metadata column that meets visualizer's requirements." >> qiime_diversity.errors   

    $qiime diversity alpha-group-significance \
        --i-alpha-diversity core-metrics-results/shannon_vector.qza \
        --m-metadata-file ${params.metadata} \
        --o-visualization core-metrics-results/shannon_group-significance.qzv 2>> qiime_diversity.errors ||\
            echo "BAC COMMENT: Could not run alpha-group-significance on shannon_vector.qza, probably due to a lack of categorical metadata column that meets visualizer's requirements." >> qiime_diversity.errors

    # This does not work if a group has only 1 member.
    # Options are to remove the group(s) with only 1 member and run the rest.
    # If not wanting default of 'permanova' then include 
    # --p-method TEXT Choices('permanova', 'anosim', 'permdisp')

    $qiime diversity beta-group-significance \
        --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
        --m-metadata-file ${params.metadata} \
        --m-metadata-column Treatment \
        --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv \
        --p-pairwise 2>> qiime_diversity.errors ||\
            echo "Could not run beta-group-significance on unweighted_unifrac_distance_matrixs.qza, probably due to a lack of categorical metadata column "Treatment" that meets visualizer's requirements." >> qiime_diversity.errors

    """
}

process get_references {
    publishDir "${baseDir}/${params.poolname}/References", mode:'copy'

    input:
    path table_art from dada2_table_for_citation
    path classifier_art from classifier_cite

    output:
    path "Software_Versions.txt"
    path "Dada2.cite"
    path "Classifier.cite"
    val 1 into completed_refs

    script:
    """
    $qiime info | awk 'NR<=5' >> Software_Versions.txt
    echo "" >> Software_Versions.txt

    $qiime cutadapt --version >> Software_Versions.txt
    echo "" >> Software_Versions.txt

    echo "Dada2 version used: " >> Software_Versions.txt
    PATH=/cluster/software/qiime2-2021.8/bin:$PATH Rscript -e "packageDescription(\\"dada2\\",fields=\\"Version\\")" | cut -d " " -f2 >> Software_Versions.txt
    echo "" >> Software_Versions.txt

    $biom --version >> Software_Versions.txt
    echo "" >> Software_Versions.txt

    echo "Classifier Info: $CLASSIFIER_INFO" >> Software_Versions.txt
    echo "" >> Software_Versions.txt

    #use this all through qiime2
    $qiime tools citations $table_art > Dada2.cite
    $qiime tools citations $classifier_art  > Classifier.cite
    """
}


process write_csv {
    publishDir "${baseDir}/${params.poolname}/", mode: 'copy'

    input:
    path rarefied_lo from rare_low_depth_csv
    path rarefied_hi from rare_hi_depth_csv
    path biom_dir from biom_make_csv

    output:
    path "results.xlsx" into results_xlsx

    script:
    """
    if [[ -f rarefied-table.high.csv ]]
    then
        tail -n +2 ${biom_dir}/annotated_table.csv > temp_ann_table.csv
        csv2xlsx_go -colsep '\\t' -overwrite --infile temp_ann_table.csv --outfile results.xlsx --append --sheet "Annotated_table" --headerlines 1
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/ASV.summary.csv --outfile results.xlsx --append --sheet "ASV_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/rep-seqs.csv --outfile results.xlsx --append --sheet "Rep_Seqs" --headerlines 0
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/Sample.summary.csv --outfile results.xlsx --append --sheet "Sample_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/stats.csv --outfile results.xlsx --append --sheet "Stats" --headerlines 2
        tail -n +2 rarefied-table.high.csv > hi_table.csv
        csv2xlsx_go -colsep '\\t' --infile hi_table.csv --outfile results.xlsx --append --sheet "Rarefied-table-high" --headerlines 1
        tail -n +2 rarefied-table.low.csv > low_table.csv
        csv2xlsx_go -colsep '\\t' --infile low_table.csv --outfile results.xlsx --append --sheet "Rarefied-table-low" --headerlines 1
    elif [[ -f rarefied-table.low.csv ]]
    then
        tail -n +2 ${biom_dir}/annotated_table.csv > temp_ann_table.csv
        csv2xlsx_go -colsep '\\t' -overwrite --infile temp_ann_table.csv --outfile results.xlsx --append --sheet "Annotated_table" --headerlines 1
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/ASV.summary.csv --outfile results.xlsx --append --sheet "ASV_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/rep-seqs.csv --outfile results.xlsx --append --sheet "Rep_Seqs" --headerlines 0
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/Sample.summary.csv --outfile results.xlsx --append --sheet "Sample_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/stats.csv --outfile results.xlsx --append --sheet "Stats" --headerlines 2
        tail -n +2 rarefied-table.low.csv > low_table.csv
        csv2xlsx_go -colsep '\\t' --infile low_table.csv --outfile results.xlsx --append --sheet "Rarefied-table-low" --headerlines 1
    else
        tail -n +2 ${biom_dir}/annotated_table.csv > temp_ann_table.csv
        csv2xlsx_go -colsep '\\t' -overwrite --infile temp_ann_table.csv --outfile results.xlsx --append --sheet "Annotated_table" --headerlines 1
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/ASV.summary.csv --outfile results.xlsx --append --sheet "ASV_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/rep-seqs.csv --outfile results.xlsx --append --sheet "Rep_Seqs" --headerlines 0
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/Sample.summary.csv --outfile results.xlsx --append --sheet "Sample_Summary"
        csv2xlsx_go -colsep '\\t' --infile ${biom_dir}/stats.csv --outfile results.xlsx --append --sheet "Stats" --headerlines 2
    fi
    """
}
