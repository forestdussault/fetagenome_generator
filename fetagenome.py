import os
import glob
import click
from subprocess import Popen, PIPE


def run_subprocess(command):
    """
    command is the command to run, as a string.
    runs a subprocess, returns stdout and stderr from the subprocess as strings.
    """
    x = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = x.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    return out, err


def kwargs_to_string(kwargs):
    """
    Given a set of kwargs, turns them into a string which can then be passed to a command.
    :param kwargs: kwargs from a function call.
    :return: outstr: A string, which is '' if no kwargs were given, and the kwargs in string format otherwise.
    """
    outstr = ''
    for arg in kwargs:
        outstr += ' {}={}'.format(arg, kwargs[arg])
    return outstr


def subsample_reads(forward_in, forward_out, num_reads, returncmd=True, reverse_in='NA', reverse_out='NA',
                    **kwargs):
    options = kwargs_to_string(kwargs)
    if os.path.isfile(forward_in.replace('_R1', '_R2')) and reverse_in == 'NA' and '_R1' in forward_in:
        reverse_in = forward_in.replace('_R1', '_R2')
        if reverse_out == 'NA':
            if '_R1' in forward_out:
                reverse_out = forward_out.replace('_R1', '_R2')
            else:
                raise ValueError('If you do not specify reverse_out, forward_out must contain _R1.\n\n')
        cmd = 'reformat.sh in1={} in2={} ' \
              'out1={} out2={} ' \
              'samplereadstarget={} {} overwrite=true'.format(forward_in, reverse_in,
                                                              forward_out, reverse_out,
                                                              str(num_reads), options)
    elif reverse_in == 'NA':
        cmd = 'reformat.sh in={} out={} samplereadstarget={} {} overwrite=true'.format(forward_in, forward_out,
                                                                                       str(num_reads), options)
    else:
        if reverse_out == 'NA':
            raise ValueError('Reverse output reads must be specified.')
        cmd = 'reformat.sh in1={} in2={} ' \
              'out1={} out2={} ' \
              'samplereadstarget={} {} overwrite=true'.format(forward_in, reverse_in,
                                                              forward_out, reverse_out,
                                                              str(num_reads), options)
    if not os.path.isfile(forward_out):
        out, err = run_subprocess(cmd)
    else:
        out = str()
        err = str()
    if returncmd:
        return out, err, cmd
    else:
        return out, err


def find_paired_reads(fastq_directory, forward_id='_R1', reverse_id='_R2'):
    """
    Looks at a directory to try to find paired fastq files. Should be able to find anything fastq.
    :param fastq_directory: Complete path to directory containing fastq files.
    :param forward_id: Identifier for forward reads. Default R1.
    :param reverse_id: Identifier for reverse reads. Default R2.
    :return: Dictionary with Sample ID as key and (R1,R2) as tupled value
    """
    pair_dict = {}
    fastq_files = glob.glob(os.path.join(fastq_directory, '*.f*q*'))
    for name in fastq_files:
        sample_id = os.path.basename(name).split(forward_id)[0]
        if forward_id in name and os.path.isfile(name.replace(forward_id, reverse_id)):
            pair_dict[sample_id] = (name, name.replace(forward_id, reverse_id))
    return pair_dict


def subsample(r1, r2, r1_out, r2_out, num_reads):
    out, err, cmd = subsample_reads(forward_in=r1, reverse_in=r2,
                                    forward_out=r1_out, reverse_out=r2_out,
                                    num_reads=num_reads)

    print('STDOUT**\n{}'.format(out))
    print('STDERR**\n{}'.format(err))
    print(cmd)
    return out, err


def normalize(r1, r2, r1_out, r2_out, target_depth=20):
    cmd = 'bbnorm.sh in={r1} in2={r2} ' \
          'out={r1_out} out2={r2_out} ' \
          'target={target_depth}'.format(r1=r1, r2=r2,
                                         r1_out=r1_out, r2_out=r2_out,
                                         target_depth=target_depth)
    out, err = run_subprocess(cmd)
    # print('STDOUT**\n{}'.format(out))
    # print('STDERR**\n{}'.format(err))
    return out, err


def subsample_folder(directory, outdir, num_reads, forward_id='R1', reverse_id='R2', ):
    """
    Normalize all of the input reads to the param provided to num_reads
    :param directory:
    :param outdir:
    :param num_reads:
    :param forward_id:
    :param reverse_id:
    :return:
    """
    fastq_dict = find_paired_reads(directory)

    subsample_dict = {}
    for key, value in fastq_dict.items():
        print('Sample ID: {}\nR1:{}\nR2:{}\n'.format(key, value[0], value[1]))
        r1, r2 = value[0], value[1]

        # Subsample read
        print('SUBSAMPLING READS...')
        r1_subsampled = os.path.join(outdir, os.path.basename(r1.replace(forward_id, forward_id + '_subsampled')))
        r2_subsampled = os.path.join(outdir, os.path.basename(r2.replace(reverse_id, reverse_id + '_subsampled')))
        subsample(r1=r1, r2=r2, r1_out=r1_subsampled, r2_out=r2_subsampled, num_reads=num_reads)

        subsample_dict[key] = (r1_subsampled, r2_subsampled)

    return subsample_dict


def count_reads(fastq_file):
    # Counts the number of reads in a .fastq.gz file
    cmd = "parallel \"gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'\" ::: "
    cmd += fastq_file
    out, err = run_subprocess(cmd)
    lines = out.splitlines()
    read_count = int(lines[0])
    return read_count


def count_reads_folder(subsampled_dict):
    read_count_dict = {}
    for key, value in subsampled_dict.items():
        read_count_dict[key] = count_reads(value[0]) + count_reads(value[1])
    return read_count_dict


def adjust_target(sample_id, subsampled_dict, read_count_dict, target_percentage, num_reads):

    # Verify read count for all samples is the same (it should be at this point)
    for key, value in read_count_dict.items():
        if int(value) != int(num_reads*2):
            print('WARNING: {} was not subsampled to the expected value of {} ({})'.format(key, num_reads, value))

    real_sampleid = str()
    for key, value in read_count_dict.items():
        if sample_id in key:
            real_sampleid = key

    total_samples = len(read_count_dict)

    # Total number of reads after removing our target sample reads
    total_reads_removed = sum(read_count_dict.values())-read_count_dict[real_sampleid]

    # Proportion of total reads that other samples occupy each
    other_percent = (1 - float(target_percentage))/(total_samples - 1)

    # Get the number of reads our target needs to be to be x% of total
    subsample_target_value = ((total_reads_removed/(other_percent*(total_samples-1))) - total_reads_removed)/2
    subsample_target_value = int(round(subsample_target_value))

    print('\nSubsampling target {} to {} reads ({}% of total of {})'.format(sample_id,
                                                                            subsample_target_value,
                                                                            float(target_percentage)*100,
                                                                            subsample_target_value+total_reads_removed))

    print(subsampled_dict[real_sampleid][0])

    r1_out = subsampled_dict[real_sampleid][0].replace('_subsampled', '_subsampled_TARGET')
    r2_out = subsampled_dict[real_sampleid][1].replace('_subsampled', '_subsampled_TARGET')

    subsample(r1=subsampled_dict[real_sampleid][0], r2=subsampled_dict[real_sampleid][1],
              r1_out=r1_out, r2_out=r2_out,
              num_reads=subsample_target_value)

    # Cleanup
    os.remove(subsampled_dict[real_sampleid][0])
    os.remove(subsampled_dict[real_sampleid][1])


def merge_fastq_dir(directory, proportion):
    cmd = 'cat *.f*q* > Merged_Fetagenome_{}.fastq.gz'.format(str(proportion))
    p = Popen(cmd, shell=True, cwd=directory)
    p.wait()


@click.command()
@click.option('-i', '--input_directory',
              type=click.Path(exists=True),
              required=True,
              help='Input directory containing all paired FASTQ files to target')
@click.option('-o', '--output_directory',
              type=click.Path(exists=True),
              required=True,
              help='Output directory to dump subsampled reads')
@click.option('-t', '--target_sample',
              required=True,
              help='ID of the target sample e.g. 2017-SEQ-0274')
@click.option('-p', '--proportion',
              required=True,
              help='Proportion of raw reads that should belong to the target sample '
                   'e.g. if set to 0.1, 10% of total reads will belong to target in the fetagenome')
@click.option('-n', '--num_reads',
              required=False,
              default=50000,
              help='Number of reads to subsample each FASTQ file to. '
                   'e.g. if set to 10, R1.fastq.gz and R2.fastq.gz will contain 10 subsampled reads each.'
                   'Default=50000.')
def fetagenome(input_directory, output_directory, proportion, target_sample, num_reads):
    # Subsample every file and dump into output directory
    print('Subsampling input directory to {} reads per read pair\n'.format(num_reads*2))
    subsample_dict = subsample_folder(directory=input_directory, outdir=output_directory, num_reads=num_reads)

    # Count number of reads for everything in the out
    read_count_dict = count_reads_folder(subsample_dict)

    # Adjust to target percentage
    adjust_target(sample_id=target_sample,
                  subsampled_dict=subsample_dict,
                  read_count_dict=read_count_dict,
                  target_percentage=proportion,
                  num_reads=num_reads)

    # Verify read counts are okay
    print('\nFINAL READ COUNTS FOR SAMPLES:')
    output_dict = find_paired_reads(output_directory)
    read_count_dict_final = count_reads_folder(output_dict)
    for key, value in read_count_dict_final.items():
        print(key + ' : ' + str(value))

    merge_fastq_dir(directory=output_directory, proportion=proportion)


if __name__ == '__main__':
    fetagenome()
