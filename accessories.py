import os
import glob
from subprocess import Popen, PIPE


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


def count_reads(fastq_file):
    # Counts the number of reads in a .fastq.gz file
    cmd = "parallel \"gunzip -c {} | wc -l | awk '{d=\$1; print d/4;}'\" ::: "
    cmd += fastq_file
    out, err = run_subprocess(cmd)
    lines = out.splitlines()
    try:
        read_count = int(lines[0])
    except IndexError:
        return None
    return read_count
