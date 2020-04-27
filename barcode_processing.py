#!/usr/bin/python3
# -*- coding: utf-8 -*-
## Author: jan@hapala.cz

"""
  Reads FASTA sequences from paired end input files and
  tries to match the left read to sequences from the reference file.
  If it fails, it tries the right read.
  Writes the statistics to the output file.

  The reference file contained three sequences (three dyes) with ambiguous nucleotides for barcoding.

  On STDOUT, it prints records for multiple files:
  -------------------------------------------------------------
  | line_start | output description                           |
  | --------------------------------------------------------- |
  | --------------------------------------------------------- |
  | ##0##      | tab-sep. barcode info for further processing |
  | ---------- | -------------------------------------------- |
  | ##1##      | similarity scores and best hit for query     |
  | ---------- | -------------------------------------------- |
  | ##2##      | failed FASTA lines (both reads)              |
  | ---------- | -------------------------------------------- |
  | ##3##      | how many 2nd reads resolved the question     |
  -------------------------------------------------------------
"""

import sys
import argparse
import pandas as pd
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


DELIMITER = "\t"
VERBOSITY = 5
MISMATCH_TOLERANCE = 0.5
FILE_STATS = "stats.csv"
FILE_BARCODES = "barcodes_clumped.csv"
FILE_FAILED_FASTA = "%s.unresolved_reads.fasta"


class BarcodeProcessor:
    """Reads DNA reads from paired FASTA files and finds best matching reference sequence.
    """

    out_dir = None
    file_stats = None
    file_barcodes = None
    file_failed_fasta = None
    ref_df = None
    delimiter = None
    logging = 0


    def __init__(self, out_dir, reference_file, logging=VERBOSITY, delimiter=DELIMITER):
        self.out_dir = out_dir
        self.file_stats = open(out_dir + "/" + FILE_STATS, "a")
        self.file_barcodes = open(out_dir + "/" + FILE_BARCODES, "a")
        self.ref_df = pd.read_csv(reference_file, sep="\t", header=0, index_col=0)
        self.logging = logging
        self.delimiter = delimiter

    def close_files(self):
        self.file_stats.close()
        self.file_barcodes.close()

    def log(self, msg_type, msg):
        """Dummy logger, prints the message to the output.

        :param msg_type: type (category or method abbreviation) of the message
        :type msg_type: str
        :param msg: message to be printed
        :type msg: str
        """
        if self.logging > 0:
            message = "## %s ## %s" % (msg_type, msg)
            print(message)

    def calculate_similarity(self, query, query_id):
        """Calculates similarity of the query and the sample.
            Applies score threshold, based on MISMATCH_TOLERANCE,
            returns only similarities ABOVE the threshold.

        :param query: DNA sequence
        :type query: str
        :param query_id: DNA sequence identificator
        :type query_id: str
        :return: similarity scores of the query sequence to all reference sequences
        :rtype: dict
        """
        similarity = { code:0 for code in self.ref_df.index.values }
        #similarity = defaultdict(0)
        print(similarity)

        for code, row in self.ref_df.iterrows(): ## iterate over ref sequences
            sequence = row.values[0] ## actual DNA sequence (code ~ seq name)
            score = 0
            
            ## reference and query are compared only if their length match
            if (len(query) == len(sequence)):
                score = sum([ 0 if sequence[i] != query[i] else 1 for i in range(len(query)) ])

            ## get threshold score for this sequence length, ignore N's
            minimum_match_score = len(sequence.replace("N", "")) * (1 - MISMATCH_TOLERANCE)
            if score < minimum_match_score: score = None
            similarity[code] = score 

        ## remove those with score None
        similarity = { code:score for code, score in similarity.items() if score is not None }

        return(similarity)


    def find_best_match(self, query, query_id, sample_name):
        """Finds the best matching sequence from the reference file
            to the query sequence

        :param query: DNA sequence
        :type query: str
        :param query_id: DNA sequence identificator
        :type query_id: str
        :param sample_name: biological sample name
        :type sample_name: str
        :return: similarity scores of the query sequence to all reference sequences
        :rtype: dict
        """
        self.log("FBM","%s" % query_id)
        similarity = self.calculate_similarity(query, query_id)

        if len(similarity) < 1:
            ## no matches for this query passing the threshold
            self.log("FBM", "RETURN NONE")
            return(None)
        else: ## success!
            ## find the key of the highest  match scoring reference sequence
            best_match_code = max(similarity.keys(), key = (lambda k: similarity[k]))
            self.log("FBM", "best_match_code %s=%s" % (best_match_code, str(similarity[best_match_code])))

            ## record the result
            self.file_barcodes.write(
                    self.delimiter.join([sample_name, query_id, best_match_code, query]) + "\n")
            self.log("FBM", "RETURN this: %s" % best_match_code)
            return(best_match_code)


    def print_reference_header(self):
        """Produces header line from the reference file
            and writes it to the output files.
            This method is run first - to initialize the output files.
        """
        keys = ["sample"] + sorted(self.ref_df.index.tolist() + ["unresolved"])
        ## print all barcodes into a header row into the stats file
        self.file_stats.write(self.delimiter.join(keys) + "\n")
        ## print header to the main output file
        self.file_barcodes.write(self.delimiter.join(["sample", "read_id", "color", "barcode"])+"\n")
    
    def print_stats(self, counts, sample_name):
        """Prints statistics

        :param counts:
        :type counts:
        :param sample_name:
        :type sample_name:
        """
        counts["unresolved"] = counts["None"]
        del counts["None"]
        codes = sorted(counts.keys())
        values = [sample_name] + [ str(counts[str(code)]) for code in codes ]
        self.log("stat", "printing... codes=%s, vals=%s" % (str(codes), str(values)))
        self.file_stats.write(self.delimiter.join(values) + "\n")
    
    def read_pair_input(self, input_file1, input_file2):
        """Reads pair DNA sequences (from two SYNCHRONIZED input FASTA files),
            finds the best matching reference sequence and returns its code.

        :param input_file1: path to input file - first reads
        :type input_file1: str
        :param input_file1: path to input file - second reads
        :type input_file1: str
        :return: generator with read pairs
        :rtype: Iterator[]
        """
        with open(input_file1, 'r') as f1, open(input_file2, 'r') as f2:

            for line_file1, line_file2 in zip(f1, f2):
                line_file1 = line_file1.strip()
                line_file2 = line_file2.strip()
                if line_file1.startswith('>'):
                    header_file1 = line_file1.replace(" ", "_")[1:]
                    header_file2 = line_file2.replace(" ", "_")[1:]
                yield(line_file1, line_file2)

    def find_pair_best_match(self, read1, header1, read2, header2, sample_name):
        """Returns the sequence code of a reference sequence best matching the query sequence,
            either its first read, or the second read if no match was found for the first.
            Returns None in case of no match.
    
        :param read1: the first read (in pair)
        :type read1: str
        :param header1: header of the first read
        :type header1: str
        :param read2: the second read (in pair)
        :type read2: str
        :param header2: header of the second read
        :type header2: str
        :param sample_name: biological sample name
        :type sample_name: str
        :return: code of the best matching reference sequence and whether the second read was matched
        :rtype: tuple (str, bool)
        """
        best_match_code = None
        is_second_read_match = False
    
        best_match_code = self.find_best_match(read1, header1, sample_name)
    
        if best_match_code is None:
            ## the 1st read in pair has not worked
            ##   -> try the 2nd one
            reverse_complement = str(Seq(read2, generic_dna).reverse_complement())
            best_match_code = self.find_best_match(reverse_complement, header2, sample_name)
            if best_match_code is not None:
                is_second_read_match = True
    
        return(best_match_code, is_second_read_match)

    def process(self, input1, input2, sample_name):
        """Reads pair DNA sequences (from two SYNCHRONIZED input FASTA files),
            finds the best matching reference sequence and returns its code.

        :param input1: path to input file - first reads
        :type input1: str
        :param input1: path to input file - second reads
        :type input1: str
        :param sample_name: biological sample name
        :type sample_name: str
        """
        failed_fasta_file_name = FILE_FAILED_FASTA % sample_name
        self.file_failed_fasta = open(self.out_dir + "/" + failed_fasta_file_name, "w")
        counts = { code:0 for code in self.ref_df.index.values }
        counts["None"] = 0
        resolved_by_pair_count = 0
        header1, header2 = "", ""

        def clean_header(header): return(header.replace(" ", "_")[1:])

        for line1, line2 in self.read_pair_input(input1, input2):

            if line1.startswith('>'):
                ## save headers
                header1, header2 = clean_header(line1), clean_header(line2)
            else:
                (best_match_code, is_second_read_match) = self.find_pair_best_match(line1, header1,
                                                                                    line2, header2, sample_name)
                if is_second_read_match: resolved_by_pair_count += 1
                if best_match_code is None:
                    self.record_non_matching(header1, line1, header2, line2)
                counts[ str(best_match_code) ] = counts[ str(best_match_code) ] + 1

        self.print_stats(counts, sample_name)
        self.file_failed_fasta.close()


    def record_non_matching(self, header1, read1, header2, read2):
        """Returns the sequence code of a reference sequence best matching the query sequence,
            either its first read, or the second read if no match was found for the first.
            Returns None in case of no match.
    
        :param read1: the first read (in pair)
        :type read1: str
        :param header1: header of the first read
        :type header1: str
        :param read2: the second read (in pair)
        :type read2: str
        :param header2: header of the second read
        :type header2: str
        """
        self.file_failed_fasta.write(">%s\n" % header1)
        self.file_failed_fasta.write("%s\n" % read1)
        self.file_failed_fasta.write(">%s\n" % header2)
        self.file_failed_fasta.write("%s\n" % read2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process FASTA sequences.')
    parser.add_argument('reffile', metavar='reffile', type=str, help='reference file')
    parser.add_argument('outdir', metavar='outdir', type=str, help='output directory')
    parser.add_argument('-1', '--infile1', metavar='infile1', type=str,
                        help='an input file1', default=None)
    parser.add_argument('-2', '--infile2', metavar='infile2',type=str,
                        help='an input file2', default=None)
    parser.add_argument('-s', '--sample_name', metavar='sample_name',type=str,
                        help='sample name', default=None)
    parser.add_argument('-v', '--verbosity', metavar='verbosity',type=int,
                        help='verbosity', default=VERBOSITY)
    parser.add_argument('-d', '--delimiter', metavar='delimiter',type=str,
                        help='delimiter', default=DELIMITER)
    args = parser.parse_args()

    processor = BarcodeProcessor(args.outdir, args.reffile,
                                 args.verbosity, args.delimiter)
    if (args.infile1 is None or args.infile2 is None or 
        args.sample_name is None): # if one of obligatory arguments is missing, print help
        processor.print_reference_header()
    else:
        processor.process(args.infile1, args.infile2, args.sample_name)
    processor.close_files()

