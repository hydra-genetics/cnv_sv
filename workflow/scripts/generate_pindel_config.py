#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

def insertSize(file):
    # Open file
    metrics = open(file)
    read_metrics = csv.reader(metrics, delimiter = "\t")

    # Loop over lines in file
    predict = []
    for row in read_metrics:
        # Consider only lines that contain 23 items
        if len(row) == 23:
            predict.append(row)
            
    # Generate dict for easy access
    parameters = dict((predict[0][n], predict[1][n]) for n in range(len(predict[0])))

    # Close file and return
    metrics.close()
    return parameters["MEAN_INSERT_SIZE"]

def writeConfigFile(output, input, insert_size, sample_id):
    with open(output[0], "wt") as out_file:
        tsv_writer = csv.writer(out_file, delimiter = "\t")
        tsv_writer.writerow([input, insert_size, sample_id])

# Call functions
writeConfigFile(snakemake.output.config, snakemake.input.bam, insertSize(snakemake.input.metrics), snakemake.wildcards.sample)
