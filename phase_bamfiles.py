#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Ryan Culligan

import os
import subprocess as sp
import pysam
import time
import gzip
import vcf

def index_samfile(sam_path):
    status = 0
    sam_index_path = sam_path + ".bai"
    if not os.path.isfile(sam_index_path):
        status = sp.call(['samtools', 'index', sam_path])
    return status


def index_all(sam_paths):
    for sam_path in sam_paths:
        index_samfile(sam_path)
    return True


def get_chrom_name(sam_path):
    return "HM763135_-_CREM_gene"
    return os.path.basename(sam_path).split(".", 1)[0].replace(" ", "_") # set chrom name



def build_paired_read_list(sam_path):
    try:
        sam_handle = pysam.AlignmentFile(sam_path, "rb")
    except ValueError:
        raise StopIteration # file doesn't have any records.  Skip.

    chrname = get_chrom_name(sam_path)
    forward_set = set()
    reverse_set = set()

    for pileupcolumn in sam_handle.pileup(chrname):
        #print (pileupcolumn.pos, pileupcolumn.n)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                if pileupread.alignment.is_read1 and (pileupread.alignment.query_name not in forward_set):
                    forward_set.add(pileupread.alignment.query_name)
                    if pileupread.alignment.query_name in reverse_set:
                        yield pileupread.alignment.query_name
                elif pileupread.alignment.is_read2 and (pileupread.alignment.query_name not in reverse_set):
                    reverse_set.add(pileupread.alignment.query_name)
                    if pileupread.alignment.query_name in forward_set:
                        yield pileupread.alignment.query_name


def build_overlap_sequences(sam_path):
    try:
        sam_handle = pysam.AlignmentFile(sam_path, "rb")
    except ValueError:
        raise StopIteration # file doesn't have any records.  Skip.

    chrname = get_chrom_name(sam_path)

    for pileupcolumn in sam_handle.pileup(chrname):

        position_set = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:

                current_base = pileupread.alignment.query_sequence[pileupread.query_position]

                if current_base not in position_set:
                    position_set[current_base] = 0
                position_set[current_base] += 1

                #print ('\tbase in read {} = {}'.format(pileupread.alignment.query_name, current_base))

        if len(position_set) > 1:
            print ("\ncoverage at base {} = {}".format(pileupcolumn.pos, pileupcolumn.n))
            print("ALLELES")
            for key in position_set:
                print("\t{}: {}".format(key, position_set[key]))


def build_paired_lists(sam_paths):
    for sam_path in sam_paths:
        print("START", sam_path)
        for i in build_paired_read_list(sam_path):
            print(i)
        print("END", sam_path)
    return True


def read_vcf(vcf_file):

    with open(vcf_file) as input_handle:
        vcf_reader =  vcf.Reader(input_handle)
        for record in vcf_reader:
            yield record


def process_vcfs_paths(vcf_paths):
    for vcf_path in vcf_paths:
        for record in read_vcf(vcf_path):
            if record.ALT != [None]:
                print(record)
        break

def fasta_diff(*fastas):
    fasta_list = []
    for fasta in fastas:
        fasta_list.append(open(fasta))

def main():
    sam_dir = "test/bams/crem"
    vcf_dir = "test/vcfs/crem"
    vcf_paths = [os.path.join(vcf_dir, i) for i in os.listdir(vcf_dir) if i.endswith(".vcf")]
    sam_paths = [os.path.join(sam_dir, i) for i in os.listdir(sam_dir) if i.endswith(".bam")]

    index_all(sam_paths)
    #build_paired_lists(sam_paths)
    build_overlap_sequences(sam_paths[-1])
    process_vcfs_paths(vcf_paths)



main()
