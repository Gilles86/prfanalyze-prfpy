#!/usr/bin/env python3
import argparse
import os
import subprocess
from glob import glob

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()


def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d" % process.returncode)


parser = argparse.ArgumentParser()

parser.add_argument('--get_config', help='Whether or not to copy the default config file.'
                    'If not specified a path to a config file is assumed to be given as second argument '
                    'and analysis will be performed.', nargs=1)

args, remainder = parser.parse_known_args()

if (args.get_config):
    copy_dir = args.get_config[0]
    run(f"cp -p default_config.yml {os.path.join(copy_dir, 'default_config.yml')}")
    run(f"chmod +w {os.path.join(copy_dir, 'default_config.yml')}")
    print(f"Config file has been copied to:\n{copy_dir}\nRestart without '--get_config' argument to start analysis.")
    exit(0)
else:
    parser = argparse.ArgumentParser()
    parser.add_argument('bids_dir', help='The directory with the input dataset '
                        'formatted according to the BIDS standard.')
    
    parser.add_argument('config_file', help='The path to the json file containing '
                        'the configuration for the fitting analysis.')
    # parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
    #                     'Multiple participant level analyses can be run independently '
    #                     '(in parallel) using the same output_dir.',
    #                     choices=['participant', 'group'])
    parser.add_argument('--output_dir', help='The directory where the output files '
                        'should be stored. If you are running group level analysis '
                        'this folder should be prepopulated with the results of the'
                        'participant level analysis.', nargs=1)
    parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                        'corresponds to sub-<participant_label> from the BIDS spec '
                        '(so it does not include "sub-"). If this parameter is not '
                        'provided all subjects should be analyzed. Multiple '
                        'participants can be specified with a space separated list.',
                        nargs="+")
    parser.add_argument('--skip_bids_validator', help='Whether or not to perform BIDS dataset validation',
                        action='store_true')
    # parser.add_argument('--debug', help='Whether to start the app in debug mode. Interactive terminal will be started.', action="store_true")
    parser.add_argument('-v', '--version', action='version',
                        version='BIDS-App example version {}'.format(__version__))


    args = parser.parse_args(remainder)

# if args.debug:
#     run('exec /bin/bash')

if not args.skip_bids_validator:
    run('bids-validator %s' % args.bids_dir)

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = {subject_dir.split(
        "-")[-1]: [] for subject_dir in subject_dirs}


# Get all sessions (not specifiable)
for subject in subjects_to_analyze.keys():
    session_dirs = glob(os.path.join(args.bids_dir, f"sub-{subject}/ses-*"))
    sessions_to_analyze = [session_dir.split("-")[-1] for session_dir in session_dirs]
    subjects_to_analyze[subject] = sessions_to_analyze

print(f"Subjects to analyze:\n{subjects_to_analyze}")

# Run analysis
for subject, sessions in subjects_to_analyze.items():
    for session in sessions:
        bold_file = os.path.join(args.bids_dir, f"sub-{subject}/ses-{session}/func/sub-{subject}_ses-{session}_task-prf_acq-normal_run-01_bold.nii.gz")
        stim_file = os.path.join(args.bids_dir, f"stimuli/sub-{subject}_ses-{session}_task-prf_apertures.nii.gz")
        # TODO the below file should be specified by the used if prfsynth has not been used
        stimjs_file = os.path.join(args.bids_dir, f"derivatives/prfsynth/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_task-prf_acq-normal_run-01_bold.json")

        if (args.output_dir):
            output_dir = os.path.join(args.output_dir[0], f"prfanalyze-prfpy/sub-{subject}/ses-{session}/")
        else:
            output_dir = os.path.join(args.bids_dir, f"derivatives/prfanalyze-prfpy/sub-{subject}/ses-{session}/")
        os.makedirs(output_dir, exist_ok=True)

        # Argument order for run_prfpy.py : (opts_file, bold_file, stim_file, stimjs_file, outdir)
        filenames_args = " ".join([args.config_file, bold_file, stim_file, stimjs_file, output_dir])
        
        run(f'python3 run_prfpy.py {filenames_args}')

