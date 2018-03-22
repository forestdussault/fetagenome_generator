import os
import click
from accessories import run_subprocess


"""
Simulating metagenomes to emulate the GRDI HiSeq output.
Read length = 150
Quality is generally very good.
Number of reads is around ~40,000,000 per sample
Looks like there is a little bit of Illumina universal adapter contamination judging by FastQC results.
"""


# NOTE: Looks like doing paired end reads makes randomreads.sh go extremely slowly. No way to multithread...
def call_randomreads(ref, output_dir, num_reads, length=150):
    r1 = os.path.join(output_dir, 'simulated_metagenomic_reads_{}.fastq.gz'.format(num_reads))
    # r2 = os.path.join(output_dir, 'simulated_metagenomic_reads_{}_R2.fastq.gz'.format(num_reads))
    cmd = 'randomreads.sh ' \
          '-Xmx40g ' \
          'ref={ref} ' \
          'out={r1} ' \
          'length={length} ' \
          'reads={num_reads} ' \
          'metagenome=t ' \
          'maxq=40 ' \
          'midq=39 ' \
          'minq=20 ' \
          'snprate=0.02 ' \
          'paired=f ' \
          'overwrite=t'.format(ref=ref, r1=r1, num_reads=num_reads, length=length)
    print(cmd)
    run_subprocess(cmd)


@click.command()
@click.option('-i', '--input_fasta',
              type=click.Path(exists=True),
              required=True,
              help='Input directory containing all paired FASTQ files to target')
@click.option('-o', '--output_dir',
              type=click.Path(exists=True),
              required=True,
              help='Directory to output simulated metagenomic reads')
@click.option('-n', '--num_reads',
              required=True,
              help='Number of reads to generate for the metagenome')
@click.option('-l', '--read_length',
              required=False,
              default=150,
              help='Length of reads to generate. Defaults to 150bp.')
def main(input_fasta, output_dir, num_reads, read_length):
    call_randomreads(ref=input_fasta, output_dir=output_dir, num_reads=num_reads, length=read_length)


if __name__ == '__main__':
    main()
