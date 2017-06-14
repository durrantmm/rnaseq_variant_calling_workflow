configfile: "config.yaml"

import os
from os.path import basename, join

WD = config['wd']

FASTQ_DIR = join(WD, config['fastq_dir'])
REFGEN_DIR = join(WD, config['refgen_dir'])
GENCODE_DIR = join(WD, config['gencode_dir'])
EXONS_BED_DIR = join(WD, config['exons_bed_dir'])
DBSNP_DIR = join(WD, config['dbsnp_dir'])
ADAR_SITES_DIR = join(WD, config['adar_sites_dir'])

STAR_GENOME_DIR = join(WD, config['star_genome_dir'])
SAM_DIR = join(WD, config['sam_dir'])
READ_GROUPS_DIR = join(WD, config['add_rg_dir'])
MARK_DUPS_DIR = join(WD, config['mark_dups_dir'])
SPLIT_READS_DIR = join(WD, config['split_reads_dir'])
RECAL_BASES_DIR = join(WD, config['recal_bases_dir'])
VARCALL_DIR = join(WD, config['varcall_dir'])
VAR_QC_DIR= join(WD, config['var_qc_dir'])
ADAR_FILT_DIR = join(WD, config['adar_filt_dir'])

WC_fastqs = glob_wildcards(join(FASTQ_DIR, '{sample}.{pair}.fq.gz'))
WC_refgens = glob_wildcards(join(REFGEN_DIR, '{genome}.fa'))
WC_gencodes = glob_wildcards(join(GENCODE_DIR, '{gencode}.gtf'))

SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']

GENOMES = set(WC_refgens.genome)
GENCODES = set(WC_gencodes.gencode)

#print(SAMPLES, PAIRS, GENOMES, GENCODES)

rule all:
    input:
        expand('{adar_filt_dir}/{{sample}}.{{genome}}.vcf'.format(adar_filt_dir=ADAR_FILT_DIR),
        sample=SAMPLES, genome=GENOMES)
    run:
        print("RNAVCW FINISHED WITH NO EXCEPTIONS!")


rule create_genome_dictionary:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR)
    run:
        command = "{picard} CreateSequenceDictionary R={input} O={output}".format(picard=config['picard_path'],
                  input=input, output=output)

        print(command)
        shell(command)


rule index_genome:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR)
    shell:
        "samtools faidx {input}"


rule bgzip_dbsnp_vcf:
    input:
        "{dbsnp_dir}/{{genome}}.vcf".format(dbsnp_dir=DBSNP_DIR)
    output:
        "{dbsnp_dir}/{{genome}}.vcf.gz".format(dbsnp_dir=DBSNP_DIR)
    shell:
        "bgzip {input}"


rule tabix_dbsnp_vcf:
    input:
        "{dbsnp_dir}/{{genome}}.vcf.gz".format(dbsnp_dir=DBSNP_DIR)
    output:
        "{dbsnp_dir}/{{genome}}.vcf.gz.tbi".format(dbsnp_dir=DBSNP_DIR)
    shell:
        "tabix -p vcf {input}"


rule make_exons_bed:
    input:
        "{gencode_dir}/{{genome}}.gtf".format(gencode_dir=GENCODE_DIR)
    output:
        "{exons_bed_dir}/{{genome}}.bed".format(exons_bed_dir=EXONS_BED_DIR)
    shell:
        "grep -P \"^[^\t]+\t[^\t]+\texon\" {input} | bedtools sort -i stdin | "
        "bedtools merge -i stdin > {output}"

rule star_genome:
    input:
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        gencode="{gencode_dir}/{{genome}}.gtf".format(gencode_dir=GENCODE_DIR)
    output:
        dir="{star_genome_dir}/{{genome}}".format(star_genome_dir=STAR_GENOME_DIR),
        sa="{star_genome_dir}/{{genome}}/SA".format(star_genome_dir=STAR_GENOME_DIR)
    shell:
        "mkdir -p {output.dir}; "
        "STAR --runThreadN 2 --runMode genomeGenerate --genomeDir {output.dir} "
        "--genomeFastaFiles {input.refgen} --sjdbGTFfile {input.gencode} --sjdbOverhang 75"

rule star_align:
    input:
        genome_dir="{star_genome_dir}/{{genome}}".format(star_genome_dir=STAR_GENOME_DIR),
        f1='{fastq_dir}/{{sample}}.R1.fq.gz'.format(fastq_dir=FASTQ_DIR),
        f2='{fastq_dir}/{{sample}}.R2.fq.gz'.format(fastq_dir=FASTQ_DIR)
    output:
        dir='{sam_dir}/{{sample}}.{{genome}}'.format(sam_dir=SAM_DIR),
        sam='{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(sam_dir=SAM_DIR)
    run:
        out_prefix = output.sam.strip('Aligned.out.sam')+'.'

        command = "mkdir -p {{output.dir}}; " \
        "STAR --readFilesIn {{input.f1}} {{input.f2}} --outFileNamePrefix {out_prefix} " \
        "--genomeDir {{input.genome_dir}} --readFilesCommand gunzip -c --runThreadN 6 " \
        "--genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 " \
        "--alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 " \
        "--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within " \
        "--outFilterType BySJout --outSAMattributes NH HI AS NM MD --sjdbScore 1 --twopassMode Basic " \
        "--twopass1readsN -1".format(out_prefix=out_prefix)

        print(command)
        shell(command)


rule add_read_groups:
    input:
        '{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.Aligned.out.sam'.format(sam_dir=SAM_DIR)
    output:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    run:
        command = "{picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id " \
                  "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample; " \
                  "echo "" > {input}".format(picard=config['picard_path'],
                  input=input, output=output)

        print(command)
        shell(command)


rule mark_duplicates:
    input:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    output:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR),
        metrics='{mark_dups_dir}/{{sample}}.{{genome}}.metrics'.format(mark_dups_dir=MARK_DUPS_DIR)
    run:
        command = "{picard} MarkDuplicates I={input} O={output_bam} CREATE_INDEX=true " \
                  "VALIDATION_STRINGENCY=SILENT M={output_metrics}".format(picard=config['picard_path'],
                  input=input, output_bam=output.bam, output_metrics=output.metrics)

        print(command)
        shell(command)


rule split_reads:
    input:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR),
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        refdict="{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR),
        refidx="{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR)
    output:
        '{split_reads_dir}/{{sample}}.{{genome}}.bam'.format(split_reads_dir=SPLIT_READS_DIR)
    run:
        command = "{gatk} -T SplitNCigarReads -I {input} -o {output} -R {refgen} " \
                  "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 " \
                  "-U ALLOW_N_CIGAR_READS".format(gatk=config['gatk_path'],
                  input=input.bam, output=output, refgen=input.refgen)

        print(command)
        shell(command)


rule recalibrate_bases:
    input:
        bam='{split_reads_dir}/{{sample}}.{{genome}}.bam'.format(split_reads_dir=SPLIT_READS_DIR),
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        dbsnp="{dbsnp_dir}/{{genome}}.vcf.gz".format(dbsnp_dir=DBSNP_DIR),
        dbsnp_idx="{dbsnp_dir}/{{genome}}.vcf.gz.tbi".format(dbsnp_dir=DBSNP_DIR),
        exons_bed="{exons_bed_dir}/{{genome}}.bed".format(exons_bed_dir=EXONS_BED_DIR)
    output:
        '{recal_bases_dir}/{{sample}}.{{genome}}.recal.table'.format(recal_bases_dir=RECAL_BASES_DIR)
    run:
        command = "{gatk} -T BaseRecalibrator -I {input} -o {output} -R {refgen} " \
                  "-knownSites {dbsnp} -L {exons_bed}".format(gatk=config['gatk_path'],
                  input=input.bam, output=output, refgen=input.refgen, dbsnp=input.dbsnp,
                  exons_bed=input.exons_bed)

        print(command)
        shell(command)


rule print_reads:
    input:
        bam='{split_reads_dir}/{{sample}}.{{genome}}.bam'.format(split_reads_dir=SPLIT_READS_DIR),
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        recal='{recal_bases_dir}/{{sample}}.{{genome}}.recal.table'.format(recal_bases_dir=RECAL_BASES_DIR)
    output:
        '{recal_bases_dir}/{{sample}}.{{genome}}.bam'.format(recal_bases_dir=RECAL_BASES_DIR)
    run:
        command = "{gatk} -T PrintReads -I {input} -o {output} -R {refgen} " \
        "-BQSR {recal}".format(gatk=config['gatk_path'],
        input=input.bam, output=output, refgen=input.refgen, recal=input.recal)

        print(command)
        shell(command)


rule haplotype_caller:
    input:
        bam='{recal_bases_dir}/{{sample}}.{{genome}}.bam'.format(recal_bases_dir=RECAL_BASES_DIR),
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        '{varcall_dir}/{{sample}}.{{genome}}.vcf'.format(varcall_dir=VARCALL_DIR)
    run:
        command = "{gatk} -T HaplotypeCaller -I {input} -o {output} -R {refgen} " \
        "-dontUseSoftClippedBases -stand_call_conf 20.0".format(gatk=config['gatk_path'],
        input=input.bam, output=output, refgen=input.refgen)

        print(command)
        shell(command)


rule variant_quality_filter:
    input:
        vcf='{varcall_dir}/{{sample}}.{{genome}}.vcf'.format(varcall_dir=VARCALL_DIR),
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        '{var_qc_dir}/{{sample}}.{{genome}}.vcf'.format(var_qc_dir=VAR_QC_DIR)
    run:
        command = "{gatk} -T VariantFiltration -V {input} -o {output} -R {refgen} " \
        "-window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD " \
        "-filter \"QD < 2.0\"".format(gatk=config['gatk_path'],
        input=input.vcf, output=output, refgen=input.refgen)

        print(command)
        shell(command)


rule make_adar_bed:
    input:
        adar_txt='{adar_sites_dir}/{{genome}}.txt'.format(adar_sites_dir=ADAR_SITES_DIR)
    output:
        adar_bed='{adar_sites_dir}/{{genome}}.bed'.format(adar_sites_dir=ADAR_SITES_DIR)
    run:
        tmp_bed = output.adar_bed+'.tmp'
        with open(input.adar_txt) as infile:
            infile.readline()
            with open(tmp_bed, 'w') as outfile:
                for line in infile:
                    line = line.strip().split()
                    chrom, start, end = line[0], int(line[1])-1, line[1]
                    outfile.write('\t'.join([chrom, str(start), end])+'\n')
        shell('bedtools sort -i {bed} > {output}; rm {bed}'.format(bed=tmp_bed, output=output.adar_bed))



rule filter_adar_sites:
    input:
        vcf='{var_qc_dir}/{{sample}}.{{genome}}.vcf'.format(var_qc_dir=VAR_QC_DIR),
        adar_bed='{adar_sites_dir}/{{genome}}.bed'.format(adar_sites_dir=ADAR_SITES_DIR)
    output:
        '{adar_filt_dir}/{{sample}}.{{genome}}.vcf'.format(adar_filt_dir=ADAR_FILT_DIR)
    shell:
        'bedtools intersect -header -v -a {input.vcf} -b {input.adar_bed} > {output}'
