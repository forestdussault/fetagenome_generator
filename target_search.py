from subprocess import Popen


def run_srst2(fastq, gene_db):
    srst2_output = fastq.replace('.fastq.gz', '_SRST2_output')
    cmd = 'srst2 ' \
          '--input_se {input_se} ' \
          '--gene_db {gene_db} ' \
          '--output {output} ' \
          '--threads {threads} ' \
          ''.format(input_se=fastq, gene_db=gene_db, output=srst2_output, threads=10)
    p = Popen(cmd, shell=True)
    p.wait()


def main():
    gene_db = '/home/forest/PycharmProjects/fetagenome_generator/amr_targets/NCBI_AMR_nt_180312_SRST2.fasta'
    fastq = '/home/forest/Documents/Projects/Fetagenome_testing/subsampled_data/subsampled_250000/merged_fetagenome.fastq.gz'

    run_srst2(fastq=fastq, gene_db=gene_db)


if __name__ == '__main__':
    main()
