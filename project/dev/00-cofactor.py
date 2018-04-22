#! /usr/bin/env python

"""
Reports the percentage of time a predicted pioneer occurs with at least one of the other target pioneers.

Helps eliminate false positives from pioneer factor prediction by finding PF cofactors, 
or TFs that bind only in the presence of at least one PF.

by Xinyu Feng, April 14 2018
"""

import sys
import pandas as pd
import numpy as np

if len(sys.argv) == 1:
	print "Usage: " + sys.argv[0] + " <pwm_ids.txt> <targets.tsv> <cell_type_1.tsv> ... <n.tsv>"
	quit()

"""
counts the number of binding and cobinding events for each target motif.
binding   := matrix[region][motif] > 0
cobinding := matrix[region][motif] > 0 and matrix[region][another motif] > 0

in:  matrix_list: list of numpy arrays 
                  (filtered matrices containing number of matches to targets)
     num_targets: total number of predicted targets
out: id2cb: column # -> # of co-binding events
     id2total: column # -> # of total binding events
"""
def cobinding(matrix_list, num_targets):

	id2cb = {}
	id2total = {}

	for i in range(num_targets):
		id2cb[i] = 0
		id2total[i] = 0

	for matrix in matrix_list:
		num_rows = matrix.shape[0]

		for ri in range(num_rows):
			row = matrix[ri]
			match_count = 0
			first_match = -1

			for mi in range(num_targets):
				if row[mi] > 0:
					id2total[mi] += 1
					if match_count == 0:
						first_match = mi 
						match_count += 1
					elif match_count == 1:
						id2cb[first_match] += 1
						id2cb[mi] += 1
						match_count += 1
					else:
						id2cb[mi] += 1
						match_count += 1

	return id2cb, id2total

## test cobinding
# test_matrix = np.array([[0, 1, 4], [3, 1, 0], [2, 0, 0], [0, 0, 1], [1, 1, 1]])
# out1, out2 = cobinding([test_matrix], 3)
# for i in range(3):
# 	print "%d: %d / %d = %f" %(i, out1[i], out2[i], float(out1[i])/out2[i])
# quit()

# load data

## targets
target_set = set()
target_list = []
id2name = {} # ex. MA0139.1	-> CTCF
num_targets = 0 # number of targets to analyze cobinding

targets_file = open(sys.argv[2])
for line in targets_file:
	tokens = line.rstrip('\n').split('\t')
	motif_id = tokens[0]
	motif_name = tokens[1]
	target_set.add(motif_id)
	target_list.append(motif_id)
	id2name[motif_id] = motif_name
	num_targets += 1
targets_file.close()

## pwms
pwm_file = open(sys.argv[1])
line_num = 0
target_cols = []
col2id = {}
for line in pwm_file:
	tokens = line.rstrip('\n').split(' ')
	motif_id = tokens[1]
	motif_name = tokens[2]
	if motif_id in target_set:
		target_cols.append(line_num)
		col2id[line_num] = motif_id
	line_num += 1 

## matrices
filtered_matrix_list = []
for i in range(3, len(sys.argv)):
	filtered_matrix = np.loadtxt(sys.argv[i], delimiter='\t', usecols=target_cols)
	filtered_matrix_list.append(filtered_matrix)

## analyze co-binding
id2cb, id2total = cobinding(filtered_matrix_list, num_targets)

## calculate rate of co-binding

for mi in range(num_targets):
	rate = float(id2cb[mi]) / id2total[mi]
	print id2name[target_list[mi]] + ': ' + str(id2cb[mi]) + ' / ' + str(id2total[mi]) + ' = ' + str(rate)
