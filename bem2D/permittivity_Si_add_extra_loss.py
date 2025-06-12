#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

take Si from BEM and add extra loss to the imaginary part of the permittivity
"""
import os

#pwd = os.getcwd()
pwd = os.path.dirname(__file__) 

#%%
# File: update_third_column.py

os.chdir(pwd)
delta = 0.1

def update_third_column(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.strip():
                continue  # Skip empty lines
            parts = line.strip().split()
            if len(parts) < 3:
                outfile.write(line)  # Leave line unchanged if less than 3 columns
                continue
            try:
                parts[2] = str(float(parts[2]) + delta)
            except ValueError:
                pass  # If third column isn't a number, leave it unchanged
            outfile.write(' '.join(parts) + '\n')

# Usage
input_filename = 'Si.txt'
output_filename = 'Si_extra_loss.txt'
update_third_column(input_filename, output_filename)
    