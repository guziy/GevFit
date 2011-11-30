Author: Huziy Oleksandr (guziy.sasha@gmail.com)

How to generate bootstrap images:
1. launch gevfit.py to calculate and save return levels
2. launch bootstrap.py to calculate fitting uncertainties
3. launch calculate_significance_from_bootstrap.py to plot
significant differences between corresponding current and future climate pairs of return levels.

There are 2 types of bootstrap that this program can perform: a) for each ensemble member separately or
b) merging all and then performing the bootstrap on merged extremes, but restricting changes to extreme values to
its corresponding member.

For b):
1. Generate return levels obtained from merged extremes: gevfit_for_all.py
2. Perform bootstrap or plot changes with significance: bootstrap_all_members_merged.py