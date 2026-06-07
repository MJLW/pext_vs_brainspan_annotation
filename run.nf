#!/usr/bin/env nextflow

params.gff = "/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/Homo_sapiens.GRCh37.87.gff3.chr.gz"
params.fasta = "/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/hs_ref_GRCh37.p5_all_contigs.fa"
params.fasta_fai = "/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/hs_ref_GRCh37.p5_all_contigs.fa.fai"
params.chr_rename = "/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/chr_rename.txt"

params.distal_tsv = "/media/data/programming/nextflow/pext_vs_bs_annotation/in/distal_regions.grch37.tsv.gz"
params.distal_hdr = "/media/data/programming/nextflow/pext_vs_bs_annotation/in/distal.hdr"

params.annotate_pext_bin = "/home/mattijn/Data/programming/c/annotate_pext/annotate_pext"
params.isoform_matrix = "/home/mattijn/Data/programming/c/annotate_pext/in/GTEx_median_tissue_expression_matrix_V8.tsv"

params.bin_and_sum_py = "/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/scripts/bin_and_sum.py"

process BCSQ {
    tag "${chr}.BCSQ"
    stageInMode "link"
    publishDir "out/bcsq", mode: "link"

    input:
    tuple val(chr), path(roulette_vcf)

    output:
    tuple val(chr), path("${chr}.bcsq.bcf")

    script:
    """
    bcftools csq --local-csq -f ${params.fasta} -g ${params.gff} ${roulette_vcf} | \
        bcftools annotate --rename-chrs ${params.chr_rename} | \
        bcftools filter \
        -i 'BCSQ[*] ~ "stop_gained" | \
            BCSQ[*] ~ "frameshift" | \
            BCSQ[*] ~ "splice_acceptor" | \
            BCSQ[*] ~ "splice_donor"' \
            -Ob -o ${chr}.bcsq.bcf
    """
}

process AnnotateDISTAL {
    tag "${chr}.DISTAL"
    stageInMode "link"
    publishDir "out/distal", mode: "link"

    input:
    tuple val(chr), path(csq_bcf)

    output:
    tuple val(chr), path("${chr}.bcsq.distal.bcf")

    script:
    """
        bcftools annotate $csq_bcf \
            -a ${params.distal_tsv} \
            -h ${params.distal_hdr} \
            -c CHROM,POS,REF,ALT,INFO/DISTAL\
            -Ob -o ${chr}.bcsq.distal.bcf
    """
}

process PEXT {
    tag "${chr}.PEXT"
    stageInMode "link"
    publishDir "out/pext", mode: "link"

    input:
    tuple val(chr), path(csq_distal_bcf)

    output:
    tuple val(chr), path("${chr}.pext.tsv")

    script:
    """
        ${params.annotate_pext_bin} ${params.isoform_matrix} $csq_distal_bcf tmp.vcf

        echo -e "CHROM\\tPOS\\tREF\\tALT\\tMR\\tGENE\\tTRANSCRIPT\\tCONSEQUENCE\\tPEXT\\tDISTAL" > ${chr}.pext.tsv
        bcftools query tmp.vcf \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MR\\t%INFO/BCSQ\\t%INFO/PEXT\\t%INFO/DISTAL\\n' | \
            awk -F'\\t' 'BEGIN{OFS="\\t"} {
                n = split(\$6, bcsq, ",");
                m = split(\$7, pext, ",");
                for (i=1; i<=n; i++) {
                    split(bcsq[i], fields, "|");
                    p = (i<=m) ? pext[i] : ".";
                    print \$1, \$2, \$3, \$4, \$5, fields[2], fields[3], fields[1], p, \$8
                }
            }' >> ${chr}.pext.tsv
    """
}

process BS {
    tag "${chr}.BS"
    stageInMode "link"
    publishDir "out/bs", mode: 'link'

    input:
    tuple val(chr), path(csq_distal_bcf)

    output:
    tuple val(chr), path("${chr}.bs.tsv")

    script:
    """
        bcftools annotate $csq_distal_bcf \
            -a {params.brainspan_tsv} \
            -h {params.brainspan_hdr} \
            -c CHROM,POS,REF,ALT,INFO/BrainSpan \
            -Ob -o tmp.bcf

        echo -e "CHROM\\tPOS\\tREF\\tALT\\tMR\\tGENE\\tTRANSCRIPT\\tCONSEQUENCE\\tBrainSpan" > ${chr}.bs.tsv
        bcftools query tmp.bcf \
            -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/MR\\t%INFO/BCSQ\\t%INFO/BrainSpan\\n' | \
            awk -F'\\t' 'BEGIN{OFS="\\t"} {
                n = split(\$6, bcsq, ",");
                m = split(\$7, brainspan, ",");
                for (i=1; i<=n; i++) {
                    split(bcsq[i], fields, "|");
                    bs = (i<=m) ? brainspan[i] : ".";
                    print \$1, \$2, \$3, \$4, \$5, fields[2], fields[3], fields[1], bs
                }
            }' >> ${chr}.bs.tsv

    """
}

process CalculateBinnedScores {
    tag "${chr}.BIN_SCORES"
    stageInMode "link"
    publishDir "out/scores/${tag}", mode: "link"

    input:
    tuple val(chr), val(tag), path(tsv)

    output:
    tuple val(chr), val(tag), path("${chr}.${tag}.tsv")

    script:
    """
        python3 ${params.bin_and_sum_py} --input $tsv --tag $tag --output ${chr}.${tag}.tsv
    """
}


workflow {
    Channel.of('chr22')
        .map { chr -> tuple(chr, file("/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/hg19_${chr}_rate_v5.2_TFBS_correction_sorted.vcf.gz")) }
        .set { roulette_ch }

    // Channel
    //    .of( (1..22).collect { "chr${it}" } + ['chrX'] )
    //    .flatten()
    //    .map { chr -> tuple(chr, file("/home/mattijn/Data/programming/nextflow/pext_vs_bs_annotation/in/hg19_${chr}_rate_v5.2_TFBS_correction_sorted.vcf.gz")) }
    //    .set { roulette_ch }

    roulette_ch | BCSQ | set { roulette_annotated }

    roulette_annotated | PEXT | map { chr, tsv -> tuple(chr, "PEXT", tsv) } | set { pext_tsvs }
    // roulette_annotated | BS | map { chr, tsv -> tuple(chr, "BrainSpan", tsv) } | set { bs_tsvs }

    pext_tsvs /*.concat(bs_tsvs)*/ | CalculateBinnedScores | set { scores }

}

