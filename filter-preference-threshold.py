#!/usr/bin/env python

import argparse
import csv
import sys

def filter_threshold(tf1_scores, tf2_scores, pref_scores, tf1_threshold, tf2_threshold):
    """
    Filters scores in a preferences bed file by returning nonzero scores only
    where the source binding probability for the TF is above that TF's threshold
    above the threshold
    :param tf1_scores: list of floats containing the scores for TF1
    :param tf2_scores: list of floats containing the scores for TF2
    :param pref_scores: list of floats containing the preference values for TF1 vs TF2
    :param tf1_threshold: NegCtrl threshold for TF1 scores
    :param tf2_threshold: NegCtrl threshold for TF2 scores
    :return: A list of preference scores, filtered by the above
    """
    output_pref_scores = list(pref_scores)
    for i in range(len(tf1_scores)):
        tf1_score, tf2_score, pref_score = tf1_scores[i], tf2_scores[i], pref_scores[i]
        if pref_score > 0.0: # favors tf_1
            # if the source TF score is below the cutoff, zero this score
            if tf1_score < tf1_threshold: output_pref_scores[i] = 0.0
        else: # favors tf_2
            if tf2_score < tf2_threshold: output_pref_scores[i] = 0.0
    return output_pref_scores


def read_scores(input, delimiter, score_index=3):
    """
    Reads the scores from a bed file into a list of floats
    :param input: Input file object
    :param delimiter: delimiter character, e.g. tab, comma, or space
    :param score_index: column index of the score value
    :return: A list of floats
    """
    scores = list()
    reader = csv.reader(input, delimiter=delimiter)
    for row in reader:
        scores.append(float(row[score_index]))
    return scores

def write_scores(scores, output, template, delimiter, score_index=3):
    """
    Writes the supplied scores to an output file
    :param scores: The scores to write into the file - list of floats
    :param output: The output file object, must be opened for writing
    :param template: The base BED file to use for writing scores. Must be opened for reading
    :param delimiter: delimiter character, e.g. tab, comma, or space
    :param score_index: column index of the score value
    """
    # Make sure we're reading from the beginning of the template output file
    template.seek(0)
    reader = csv.reader(template, delimiter=delimiter)
    writer = csv.writer(output, delimiter=delimiter)
    for i, row in enumerate(reader):
        row[score_index] = scores[i]
        writer.writerow(row)

def main(args):
    tf1_scores = read_scores(args.tf1_bedfile, args.delimiter)
    tf2_scores = read_scores(args.tf2_bedfile, args.delimiter)
    pref_scores = read_scores(args.pref_bedfile, args.delimiter)
    lengths_equal = len(tf1_scores) == len(tf2_scores) == len(pref_scores)
    if not lengths_equal:
        raise RuntimeError('Input files must all be same length')
    output_pref_scores = filter_threshold(tf1_scores, tf2_scores, pref_scores, args.tf1_threshold, args.tf2_threshold)
    args.tf1_bedfile.seek(0)
    write_scores(output_pref_scores, sys.stdout, args.tf1_bedfile, args.delimiter)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Filter preferences scores based on input prediction values and NegCtrl thresholds')
    parser.add_argument('tf1_bedfile', type=argparse.FileType('rb'))
    parser.add_argument('tf2_bedfile', type=argparse.FileType('rb'))
    parser.add_argument('pref_bedfile', type=argparse.FileType('rb'))
    parser.add_argument('tf1_threshold', type=float)
    parser.add_argument('tf2_threshold', type=float)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--spaces', action='store_const', const=b' ', dest='delimiter')
    group.add_argument('--tabs', action='store_const', const=b'\t', dest='delimiter')
    group.add_argument('--commas', action='store_const', const=b',', dest='delimiter')
    parser.add_argument_group()
    args = parser.parse_args()
    main(args)

