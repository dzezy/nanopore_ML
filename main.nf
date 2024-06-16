#!/usr/bin/env nextflow

// merge multiple CRAM/BAM files, truncate output filename to HG00X_merged.bam
process merge_crams {
    conda "$projectDir/envs/samtools.yaml"

    cache 'lenient'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_merged.bam")

    cpus params.threads

    """
    samtools merge -f -@ ${task.cpus} -o ${sample_id}_merged.bam ${bam}
    """
}

// convert BAM to FASTQ format for remapping
process bamtofastq {
    conda "$projectDir/envs/samtools.yaml"

    cache 'lenient'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq")

    cpus params.threads

    """
    samtools fastq -@ ${task.cpus} ${bam} > ${sample_id}_merged.fastq
    """
}

// do QC on fastq reads with fastqc
process fastqc {
    publishDir "$projectDir/output/fastqc", mode: 'copy'

    conda "$projectDir/envs/fastqc.yaml"
    cache 'lenient'

    input:
    tuple val(sample_id), path(fastq)

    output:
    path "${sample_id}_merged_fastqc.html"
    path "${sample_id}_merged_fastqc.zip"

    cpus params.threads

    """
    fastqc -t ${task.cpus} --memory 4000 ${fastq}
    """
}

// map nanopore fastq reads with minimap2
process minimap2 {
    conda "$projectDir/envs/minimap2.yaml"

    cache 'lenient'

    publishDir "$projectDir/output/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}_merged_realigned.bam")

    cpus params.threads

    """
    minimap2 -ax map-ont -t ${task.cpus} ${params.reference} ${fastq} | samtools view -hb > ${sample_id}_merged_realigned.bam
    """
}


// generate mapping statistics with flagstat
process flagstat {
    conda "$projectDir/envs/samtools.yaml"

    cache 'lenient'

    publishDir "$projectDir/output/flagstat", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("${bam.baseName}.flagstat.txt")

    cpus params.threads

    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${bam.baseName}.flagstat.txt
    """
}



// sort and index CRAM/BAM file
process sort_index {
    conda "$projectDir/envs/samtools.yaml"

    cache 'lenient'

    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${bam.baseName}.sorted.bam"), path("${bam.baseName}.sorted.bam.bai")

    cpus params.threads

    """
    samtools sort -o ${bam.baseName}.sorted.bam ${bam}
    samtools index -@ ${task.cpus} ${bam.baseName}.sorted.bam
    """
}

// variant calling on CRAM/BAM file with bcftools mpileup and call
process variant_call {
    conda "$projectDir/envs/lofreq.yaml"

    cache 'lenient'

    publishDir "$projectDir/output/unfiltered_vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(bed)

    output:
    tuple val(sample_id), path("${bam.simpleName}_highconf.vcf")

    cpus params.threads
    
    """
    lofreq call-parallel --pp-threads ${task.cpus} --no-default-filter -B -q 3 -Q 3 -a 1 -b 1 -f ${params.reference} -l ${bed} -o ${bam.simpleName}_highconf.vcf ${bam}
    """
}

// filter out variants that have less than a given minimum alt coverage 
process variant_filter {
    conda "$projectDir/envs/bcftools.yaml"

    cache 'lenient'

    publishDir "$projectDir/output/filtered_vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${vcf.simpleName}_filtered_altcov_${params.min_alt_cov}.vcf")

    """
    bcftools filter -i '(DP4[3]+DP4[4]) >= ${params.min_alt_cov}' ${vcf} -o "${vcf.simpleName}_filtered_altcov_${params.min_alt_cov}.vcf"
    """
}

// extract X bp flanking sequences for variant sites in vcf file
process extract_seq_context {
    conda "$projectDir/envs/bedtools.yaml"

    cache 'lenient'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${vcf.baseName}_flanking_seq.tsv")

    """
    vcf2bed --do-not-split < ${vcf} | awk '{print \$1"\t"\$2"\t"\$3}' | \
    bedtools slop -i - -g ${params.fai} -b ${params.flanking_seq_length} | \
    bedtools getfasta -tab -fi ${params.reference} -bed - -fo ${vcf.baseName}_flanking_seq.tsv
    """
}

// parse sequence context tsv files in preparation for dataframe reading
process parse_seq_context {
    cache 'lenient'

    publishDir "$projectDir/output/seq_context", mode: 'copy'

    input:
    tuple val(sample_id), path(tsv)

    output:
    tuple val(sample_id), path("${tsv.baseName}_parsed.tsv")

    shell:
    '''
    awk -F'\t' '{
        split($1, chr_parts, /[:-]/)

        $1 = chr_parts[1]
        $3 = $2
        $2 = chr_parts[3] - 100

        # Print the modified line
        print
    }' !{tsv} > !{tsv.baseName}_parsed.tsv
    '''
}

// concatenate tsv sequences and vcf together for csv reading
process concat_seqs {
conda "$projectDir/envs/bedtools.yaml"

    cache 'lenient'

    publishDir "$projectDir/output/final_tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tsv)

    output:
    path("${vcf.baseName}_with_seq.tsv")

    """
    grep -v "#" ${vcf} | sort -k 1,1 -k2,2n > ${sample_id}_intermediate.vcf
    grep -v "#" ${tsv} | sort -k 1,1 -k2,2n > ${sample_id}_intermediate.tsv

    paste -d'\t' ${sample_id}_intermediate.vcf <(cut -d' ' -f3 ${sample_id}_intermediate.tsv) > ${vcf.baseName}_with_seq.tsv
    """
}


// define the workflow steps
workflow {

    pairs = Channel.fromFilePairs(params.cram)
    pairs
        .map{ it -> [it[0].substring(0,5), it[1][0]] }
        .groupTuple()
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .set {grouped_samples}
        
    bams = merge_crams(grouped_samples) | bamtofastq | minimap2 | sort_index
    mapping_stat = flagstat(minimap2.out)

    beds = Channel.of(['HG001', params.bed1], ['HG002', params.bed2], ['HG003', params.bed3], ['HG004', params.bed4])
    
    bams
        .join(beds)
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .set {bam_with_bed}

    vcfs = variant_call(bam_with_bed) | variant_filter
    giabs = Channel.of(['HG001_GIAB', params.giab1], ['HG002_GIAB', params.giab2], ['HG003_GIAB', params.giab3], ['HG004_GIAB', params.giab4])
    mixed_vcfs = vcfs.concat(giabs)
    context = extract_seq_context(mixed_vcfs) | parse_seq_context
    
    result_tsv = mixed_vcfs.join(context) | concat_seqs
    result_tsv.view()
}
