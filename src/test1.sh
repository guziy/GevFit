#!/bin/bash
##$ -j y
#$ -cwd
#$ -e /home/huziy/skynet3_rech1/error
#$ -o /home/huziy/skynet3_rech1/log
#$ -S /bin/sh
#$ -M guziy.sasha@gmail.com
#$ -m b
#$ -q q_skynet3
##ls -l>/home/nadjet/list
#$ -pe shm 20
##Use my current environment variables
#$ -V
which python
python bootstrap_all_members_merged.py 
