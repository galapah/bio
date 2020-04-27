#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
## Parser of FASTA files, splits the input into two FASTA files
##  based on a gene code list (selected sequences)
##
## Author: jan@hapala.cz
##
## USAGE EXAMPLE
##
# ./separate_sequences.py ../annot/TAIR10_downstream_500_20101028.fa ../annot/non_protein_codes.txt ../annot/TAIR10_downstream_500_20101028_NON_PROTEIN.fa ../annot/TAIR10_downstream_500_20101028_PROTEIN.fa


import os, re, sys, doctest, time


class FASTAParserSeparator():
    """
    """
    results = []
    ifilename = None
    ofilename_pos = None
    ofilename_neg = None
    codes = None


    def __init__(self, ifilename, codes_filename, 
                       ofilename_pos, ofilename_neg):
        """
        :param ifilename: input FASTA file
        :type ifilename: str
        :param codes_filename: file with list of gene codes
        :type codes_filename: str
        :param ofilename_pos: output FASTA file with sequences from the list
        :type ofilename_pos: str
        :param ofilename_neg: output FASTA file with sequences NOT in the list
        :type ofilename_neg: str
        """
        self.ifilename = ifilename
        self.ofilename_pos = ofilename_pos
        self.ofilename_neg = ofilename_neg
        self.codes = self.load_codes(codes_filename)

    def load_codes(self, codes_file):
        """Loads gene codes from file and returns
            codes without duplicates.

        :param codes_file: path to the list of codes
        :type codes_file: str
        :return: gene codes without duplicates
        :rtype: set
        """
        with open(codes_file, "r") as cfile:
            codes = set([ line.strip().lower() for line in cfile ])
        if "" in codes: codes.remove("")
        return(codes)

    def next_sequence(self):
        """Generator
            Reads FASTA input and serves it by one sequence,
            triplet (code, header, sequence)

        :return: triplet (code, header, sequence)
        :rtype: Iterator[]
        """
        code, sequence = "", ""
        with open(self.ifilename, "r") as fr:
            for line in fr:
                line = line.strip()
                if line.startswith(">"): # is header line
                    if code != "":
                        # new sequence encountered, serve the previous one
                        yield(code, header, sequence)
                    header, code, sequence = line, _extract_code(line), sequence
                else:
                    sequence += line
            # serve the last sequence
        yield(code, header, sequence)

    def parse(self):
        """Loops through the input FASTA file,
            checks their gene codes - if they are on the list.
            Based on that, splits it into two output FASTA files.
        """
        with open(self.ofilename_pos, "w") as ofile_pos,\
             open(self.ofilename_neg, "w") as ofile_neg:
    
            for (count, (code, header, sequence)) in enumerate(
                                            self.next_sequence()):
                # pick the right output file
                if code in self.codes: ofile = ofile_pos
                else: ofile = ofile_neg
                ofile.write(f"{header}\n{sequence}\n")


def _extract_code(header):
    """Extracts gene code from FASTA header line

    :param header: FASTA sequence header
    :type header: str
    :return: gene code of the FASTA sequence
    :rtype: str

    >>> _extract_code(">BRCA2 gene involved in breast cancer")
    'brca2'
    >>> _extract_code("> BRCA2 gene involved in breast cancer")
    ''
    """
    code = ""
    matcher = re.compile(">([^\s]*).*")
    m = matcher.match(header)
    try:
        code = m.group(1)
    except AttributeError as e:
        raise ValueError(f"No code in header: {header}") from e
    return (code.lower())


if __name__ == "__main__":
    doctest.testmod()
    args = sys.argv[1:]
    if len(args) < 4:
        raise Exception("Not enough arguments, dear!")

    #path = args[0]
    ifilename, codes_filename = args[0], args[1]
    ofilename_pos, ofilename_neg = args[2], args[3]
    parser = FASTAParserSeparator(ifilename, codes_filename,
                                  ofilename_pos, ofilename_neg)
    parser.parse()

