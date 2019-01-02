#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Written by Shikai Jin on 2018-Jun-18, latest modified on 2018-Jun-28

import argparse
import subprocess


def submit_jobs(index, jobnum):
    slurm = 'job' + str(index) + '.pbs'
    if index == 1:
        cmd = 'sbatch ' + slurm
    else:
        cmd = 'sbatch --depend=afterany:' + jobnum + ' ' + slurm
    print ("The code has submitted job with command: %s" % (cmd))

    (status, jobnum) = subprocess.getstatusoutput(cmd)

    if (status == 0):
        jobnum = jobnum.split()[3]
        print ("Job number is %s" % (jobnum))
    else:
        print ("Error! The status code is %s" % (status))
        print ("The reason is:\n%s" % (jobnum))
        exit()
    return jobnum


def main():
    parser = argparse.ArgumentParser(
        description="This script submits jobs on Slurm, you can choose single job or multiple interdependece jobs")
    parser.add_argument("rounds", help="How many jobs you want", type=int)
    args = parser.parse_args()
    rounds = args.rounds

    jobnum = 0
    for index in range(1, rounds + 1):
        if rounds == 1:
            (status, jobnum) = subprocess.getstatusoutput('sbatch job.pbs')
            if (status == 0):
                jobnum = jobnum.split()[3]
                print ("Job number is %s" % (jobnum))
            else:
                print ("Error! The status code is %s" % (status))
                print ("The reason is:\n%s" % (jobnum))
                exit()
        else:
            jobnum = submit_jobs(index, jobnum)


if __name__ == '__main__':
    main()
