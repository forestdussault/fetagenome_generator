import os
import glob
import click
import pandas as pd
from subprocess import Popen
from fetagenome import find_paired_reads


def run_srst2(sampleid, r1, r2, gene_db, output_directory, threads=10):
    """
    Calls SRST2 on paired-end reads
    :param sampleid: Base sample ID for reads
    :param r1: Forward read (.fastq.gz)
    :param r2: Reverse read (.fastq.gz)
    :param gene_db: FASTA file containing AMR targets
    :param output_directory: Destination of output files
    :param threads: Number of CPUs to dedicate to task
    """
    print('\nRunning SRST2 on {}'.format(sampleid))
    srst2_output = os.path.join(output_directory, sampleid + '_SRST2')
    cmd = 'srst2 ' \
          '--input_pe {r1} {r2} ' \
          '--gene_db {gene_db} ' \
          '--output {output} ' \
          '--threads {threads} ' \
          ''.format(r1=r1, r2=r2, gene_db=gene_db, output=srst2_output, threads=threads)
    p = Popen(cmd, shell=True)
    p.wait()

    # Remove extra junk
    txt_files = glob.glob(os.path.join(output_directory, '*.txt'))
    txt_files = [x for x in txt_files if 'fullgenes' in x]
    for file in os.listdir(output_directory):
        file = os.path.join(output_directory, file)
        if file not in txt_files:
            os.remove(file)


def combine_results(srst2_result_directory):
    """
    Combines all SRST2 fullgenes file output into one tab delimited file
    :param output_directory: Directory containing SRST2 output
    """
    output_files = glob.glob(os.path.join(srst2_result_directory, '*fullgenes*'))

    df_list = []
    for file in output_files:
        tmp_df = pd.read_csv(file, delimiter='\t')
        df_list.append(tmp_df)

    df = pd.concat(df_list, ignore_index=True)
    df.to_csv(os.path.join(srst2_result_directory, 'SRST2_ResFinder_Combined_Output.tsv'), sep='\t')


@click.command()
@click.option('-i', '--input_directory',
              type=click.Path(exists=True),
              required=True,
              help='Input directory containing all paired FASTQ files to target')
@click.option('-o', '--output_directory',
              type=click.Path(exists=True),
              required=True,
              help='Output directory to dump SRST2 results')
@click.option('-g', '--gene_db',
              type=click.Path(exists=True),
              required=True,
              help='Gene targets to search against')
def main(input_directory, output_directory, gene_db):
    # Get pairs
    pair_dict = find_paired_reads(input_directory)

    # Run on every sample
    for sampleid, reads in pair_dict.items():
        run_srst2(sampleid=sampleid,
                  r1=reads[0],
                  r2=reads[1],
                  gene_db=gene_db,
                  output_directory=output_directory)

    # Create final output file
    combine_results(srst2_result_directory=output_directory)


if __name__ == '__main__':
    main()
