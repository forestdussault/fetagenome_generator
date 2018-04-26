import os
import glob
import click
from accessories import run_subprocess


# NOTE: Looks like doing paired end reads makes randomreads.sh go extremely slowly. No way to multithread...
def call_randomreads(ref, output_dir, num_reads, length=150):
    r1 = os.path.join(output_dir, 'Simulated_Metagenome_{}reads.fastq.gz'.format(num_reads))
    # r2 = os.path.join(output_dir, 'simulated_metagenomic_reads_{}_R2.fastq.gz'.format(num_reads))
    cmd = 'randomreads.sh ' \
          '-Xmx40g ' \
          'ref={ref} ' \
          'out={r1} ' \
          'length={length} ' \
          'reads={num_reads} ' \
          'metagenome=t ' \
          'maxq=40 ' \
          'midq=30 ' \
          'minq=20 ' \
          'snprate=0.02 ' \
          'paired=f ' \
          'overwrite=t'.format(ref=ref, r1=r1, num_reads=num_reads, length=length)
    print(cmd)
    out, err = run_subprocess(cmd)
    print(out)
    print(err)
    return r1


def fuse_contigs():
    pass


def concatenate_fasta_list(input_dir):
    fasta_list = glob.glob(os.path.join(input_dir, '*.f*a'))
    fasta_string = str()
    output_fasta = os.path.join(os.path.join(input_dir, 'ConcatenatedAssemblies_temp.fasta'))
    for fasta in fasta_list:
        fasta_string += fasta + ' '
    cmd = 'cat {} > {}'.format(fasta_string, output_fasta)
    run_subprocess(cmd)
    return output_fasta


@click.command()
@click.option('-i', '--input_dir',
              type=click.Path(exists=True),
              required=True,
              help='Input directory containing all assemblies to use as a source for your simulated metagenome')
@click.option('-o', '--output_dir',
              type=click.Path(),
              required=True,
              help='Directory to output simulated metagenomic reads')
@click.option('-n', '--num_reads',
              required=True,
              help='Number of reads to generate for the metagenome')
@click.option('-l', '--read_length',
              required=False,
              default=150,
              help='Length of reads to generate. Defaults to 150bp.')
def main(input_dir, output_dir, num_reads, read_length):
    # setup
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # run program
    concatenated_fasta = concatenate_fasta_list(input_dir)
    metagenome = call_randomreads(ref=concatenated_fasta, output_dir=output_dir,
                                  num_reads=num_reads, length=read_length)
    print('Simulated metagenome path: {}'.format(metagenome))

    # cleanup
    os.remove(concatenated_fasta)


if __name__ == '__main__':
    main()
