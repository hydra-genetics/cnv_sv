#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@igp.uu.se"
__license__ = "GPL-3"


import csv
import logging


def insertSize(file):
    # Open file
    metrics = open(file)
    read_metrics = csv.reader(metrics, delimiter="\t")

    # Loop over lines in file
    predict = []
    for row in read_metrics:
        # Consider only lines that contain 23 items
        if len(row) == 23:
            predict.append(row)
    metrics.close()

    # Generate dict for easy access
    try:
        parameters = dict(
            (predict[0][n], predict[1][n])
            for n in range(len(predict[0]))
        )
    except IndexError:
        logging.error("Could not detect insert size in %s", file)
        raise

    # Return insert size
    try:
        parameters["MEAN_INSERT_SIZE"]
    except KeyError:
        logging.error("Could not detect insert size in %s", file)
        raise
    else:
        return parameters["MEAN_INSERT_SIZE"]


def writeConfigFile(output, input, insert_size, sample_id):
    with open(output, "wt") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        tsv_writer.writerow([input, insert_size, sample_id])


# Call functions
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, filename=snakemake.log[0])
    logging.info("Read %s", snakemake.input.metrics)
    insert_size = insertSize(snakemake.input.metrics)
    logging.info("Insert size determined to be %s", insert_size)
    writeConfigFile(
        snakemake.output.config,
        snakemake.input.bam,
        insert_size,
        snakemake.wildcards.sample,
    )
    logging.info(
        "Successfully written config file %s",
        snakemake.output.config,
    )
