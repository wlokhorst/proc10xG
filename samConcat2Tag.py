#!/usr/bin/env python
'''
Copyright 2017 Matt Settles
Created June 8, 2017

Revised by Wouter Lokhorst, 5FEB2019

bwa mem -C option concatenates the fasta/fastq
CTACATTGTCAAGGGT:E00558:34:HGCJ3ALXX:1:1101:2108:1731   99      000000F 9225\
71  60      127M    =       922961  517     ACTCGGGGAGGTGTTAGCTGCTGCCTCACACA\
TTGGGTTTATAGGCTGAATCTTGTTCTCTTTAGGCTTCCAGAGTTTTCTCAGTTACTATTTCTCCTGTCACATACT\
CGCTGCTTCTTCTGTCATA JJJJJJ<JJF<7A7FJJJJJJ<JJJAJAJJJFJFFFJ----AJJFJ---7---<FJ\
J<JF<7FFFJJJFJJAJF-AAFFFFF-AFJF7FF<A--FJJJAF)-7-77<<7--)7)<<--77A7-<--< NM:i\
:3  MD:Z:74T34A3T13 AS:i:112        XS:i:19 1:N:0:GOOD:CCGATTAA:CTACATTGTCAA\
GGGT:<AAFFJJFJJFJJJJJ:CCAGTGA:J<FFFJJ

This pulls it out, 9 columns and produces new 10x tags in the bam then write\
s to out
'''
import sys
import os
import argparse
from itertools import islice
import multiprocessing as mp

lock = mp.Lock()


def make_output_file(out_name):
    """Creates the output file.
    """
    with open(out_name, "w") as outfile:
        pass


def write_line(out_name, line):
    """Writes a single line to either stdout or the output file.
    """
    lock.acquire()
    # Checks if out_name is a string.
    if isinstance(out_name, str):
        out = open(out_name, "a")
        out.write(line)
        out.close()
    else:
        sys.stdout.write(line)
    lock.release()


def extract_tag(line, out_name):
    """Extracts the GEM tag for a single alignment line.

    Input:
    - line: string (the alignment line)
    - out_name: either a string (full path of the output file) or stdout
    Output:
    - An alignment file without the GEM tags
    """
    line2 = line.strip().split()
    # Comment/header lines start with @
    if line[0] != "@" and len(line2) > 2:
        # Handle TAG:
        # get the final concatenated tag
        tag = line2[-1]
        if (tag[0:6] in ['1:N:0:', '2:N:0:']):
            tsplit = tag.split(":", 4)
            tsplit2 = tsplit[4].split("_")
            if len(tsplit2) != 5:
                sys.stderr.write("SAMCONCAT\tERROR\tsam file has \
concatenated info, but its the wrong size")
                sys.exit(1)
            # fixed barcode
            line2 = line2[0:-1]
            line2.extend(["ST:Z:" + tsplit2[0],
                          "BX:Z:" + line2[0].split(":")[0] + "-1",
                          "BC:Z:" + tsplit[3],
                          "QT:Z:" + '!' * len(tsplit[3]),
                          "RX:Z:" + tsplit2[1],
                          "QX:Z:" + tsplit2[2],
                          "TR:Z:" + tsplit2[3],
                          "TQ:Z:" + tsplit2[4]])
            write_line(out_name, '\t'.join(line2) + '\n')
        else:   # Does not contain a concatenated tag as expected by bwa mem
            write_line(out_name, line)
    else:  # Its the header lines, so just put back on the stream/file
        write_line(out_name, line)


class ArgumentParserError(Exception):
    """Creates a new kind of error.
    """
    pass


class ThrowingArgumentParser(argparse.ArgumentParser):
    """Redefines the argparse error, to be able to catch it with try/except.
    """
    def error(self, message):
        raise ArgumentParserError(message)


def handle_args():
    """Handles arguments both in the command line and in IDLE.
    
    Output:
    Tuple, consisting of:
    - string (input filename or stdin)
    - string (output filename or stdout)
    - integer (number of CPUs)
    """
    version_num = "0.0.2"
    # Tries to execute the script with command line arguments.
    try:
        # Creates an instance of argparse.
        argparser = ThrowingArgumentParser(prog=sys.argv[0],
            description='samConcat2Tag, processes bwa mem sam format where \
the read comment has been appended to the mapping line following process_10\
xReads.py', epilog='For questions or comments, please contact Matt Settles \
<settles@ucdavis.edu>\n%(prog)s version: ' + version_num, add_help=True)
    except ArgumentParserError:
        print("Please run this script on the command line, with the \
correct arguments. Type -h for help.\n")
        sys.exit()
    else:
        # Adds the positional arguments.
        argparser.add_argument('inputfile', metavar='inputsam', type=str,
            nargs='?', help='Sam file to process [default: %(default)s]',
            default="stdin")
        # Adds the optional arguments.
        argparser.add_argument('--version', action='version',
            version="%(prog)s version: " + version_num)
        # TODO: ADD parameter for sample ID
        argparser.add_argument('-o', '--output_base',
            help="Directory + prefix to output, [default: %(default)s]",
            action="store", type=str, dest="output_base", default="stdout")
        argparser.add_argument("-@", "--cpus",
            help="The number of CPUs to use.", type=int, default=1)
        # Parses the arguments given in the shell.
        args = argparser.parse_args()
        inp = args.inputfile
        outb = args.output_base
        cpus = args.cpus
    return inp, outb, cpus


if __name__ == "__main__":
    inp, outb, cpus = handle_args()
    if outb == "stdout":
        out = FALSE
    else:
        out = outb + ".sam"
        make_output_file(out)
    if inp == 'stdin':
        # Reading from stdin.
        insam = sys.stdin
    else:
        if not os.path.exists(inp):
            sys.exit("Error, can't find input file %s" % inp)
        insam = open(inp, 'r')
    # Maximum number of concurrent processes is the given number of CPUs.
    P = mp.Pool(cpus)
    # Read the file line by line, without loading it all into memory.
    while True:
        chunk = list(islice(insam, 1))
        if not chunk:
            break
        line = chunk[0]
        P.apply_async(extract_tag,args=(line, out,))
    P.close()
    P.join()
    # Checks if insam is a string.
    if isinstance(insam, str):
        insam.close()
