#! /usr/bin/env python

import click
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# Usage something like:
# makeShortContigs.py -l 200000 -c 1,2,3,X -r human_g1k_v37_decoy.fasta

@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--length',      '-l', type=int, help='length of the contigs to export', required=True)
@click.option('--contigs',     '-c', type=str, help='ids of contigs in the fasta file', required=True)
@click.option('--reference',   '-r', type=str, help='source reference file', required=True)
def exportShortContigs(length,contigs,reference):
    # this is the main processing routine
	contigsToPrint = contigs.split(",")
	for seq_record in SeqIO.parse(reference,"fasta"):
		# seq_record.id is something like
		# >chr1 dna:chromosome chromosome:GRCh37:1:1:249250621:1
		# but we want to have the "chr1" only, so have to split and replace
		shortID = seq_record.id.split()[0].replace(">","")
		if shortID in contigs:
			newSeq = seq_record.seq[0:length]
			sys.stdout.write( SeqRecord(newSeq, id=seq_record.id, description=seq_record.description).format("fasta") )

if __name__ == "__main__":
    exportShortContigs()

