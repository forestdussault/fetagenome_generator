import os
import subprocess

"""
This file will allow you to run the fetagenome.py and srst2_wrapper.py scripts for each proportion specified
multiplied by the number of replicates
"""


def main():
    # Where the folder structure for the analysis will be created
    base_folder = '/home/forest/Documents/Projects/Fetagenome/test'

    # Where all of your raw read files are located
    data_folder = '/home/forest/Documents/Projects/Fetagenome/raw_data'

    # Your AMR database location
    srst2_targets = '/home/forest/PycharmProjects/fetagenome_generator/amr_targets/NCBI_AMR_nt_180312_SRST2.fasta'

    # The target proportions
    proportions = [
        '0.1',
        '0.08',
        '0.06',
        '0.04',
        '0.02',
        '0.01',
        '0.005'
    ]

    # Number of replicates you'd like to do
    num_replicates = 10

    # Create base analysis folder if it doesn't already exist
    if not os.path.isdir(base_folder):
        os.mkdir(base_folder)

    # Run scripts for each proportion*number_replicates
    for proportion in proportions:
        # make directory
        try:
            os.mkdir(os.path.join(base_folder, proportion))
        except OSError:
            pass

        # 10 replicates
        replicate_paths = []
        for x in range(num_replicates):
            p = os.path.join(base_folder, proportion, 'rep' + str(x+1))

            # make directory
            try:
                os.mkdir(p)
            except OSError:
                pass
            replicate_paths.append(p)

        for output_path in replicate_paths:
            # Generate fetagenome with the following parameters.
            # Adjust -n parameter to the minimum number of reads in your FASTQ raw dataset
            # Set -t to the target SeqID you're adjusting the proportion of
            cmd = 'python fetagenome.py ' \
                  '-i {input} ' \
                  '-o {output} ' \
                  '-n 300000 ' \
                  '-t 2016-SEQ-1009 ' \
                  '-p {proportion}'.format(input=data_folder, output=output_path, proportion=proportion)
            p = subprocess.Popen(cmd, shell=True)
            p.wait()

            # Run SRST2
            cmd_srst2 = 'python srst2_wrapper.py ' \
                        '-i {output_path} ' \
                        '-o {output_path} ' \
                        '-g {srst2_targets} ' \
                        '-se'.format(srst2_targets=srst2_targets, output_path=output_path)
            p = subprocess.Popen(cmd_srst2, shell=True)
            p.wait()

    print('FINISHED RUNNING')


if __name__ == '__main__':
    main()
