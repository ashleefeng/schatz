#! /usr/bin/env python

import numpy as np 
import pandas as pd 
import sys

filenames = sys.argv[1:]

# results = pd.DataFrame(columns=["q1.1", "q1.2", "q1.3a", "q1.3b", "q1.4a", "q1.4b", "q1.5"])

index = ["Total genome size", "Number of chromosomes", "Largest chromosome size", "Largest chromosome name", \
"Smallest chromosome size", "Smallest chromosome name", "Mean chromosome length"]

results = pd.DataFrame(index = index)


for filename in filenames:
	print "Parsing", filename
	tokens = filename.split('/')[-1].split('.')
	species = tokens[0]
	df = pd.read_csv(filename, sep='\t', header=None, names=["chr", "size"])
	total_size = df["size"].sum(0)
	nchrom = len(df)
	max_chrom_size = max(df["size"])
	max_chrom_name = str(df["chr"][np.argmax(df["size"])])
	min_chrom_size = min(df["size"])
	min_chrom_name = str(df["chr"][np.argmin(df["size"])])
	mean_chrom_len = df["size"].mean(0)

	results[species] = [total_size, nchrom, max_chrom_size, max_chrom_name, \
	                    min_chrom_size, min_chrom_name, mean_chrom_len]


results.to_csv("Q1.csv")