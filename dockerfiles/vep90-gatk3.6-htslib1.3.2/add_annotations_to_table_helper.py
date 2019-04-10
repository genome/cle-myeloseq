#!/usr/bin/env python

import sys
import os
import re
from cyvcf2 import VCF
import tempfile
import csv
import binascii

def parse_csq_header(vcf_file):
    for header in vcf_file.header_iter():
        info = header.info(extra=True)
        if b'ID' in info.keys() and info[b'ID'] == b'CSQ':
            format_pattern = re.compile('Format: (.*)"')
            match = format_pattern.search(info[b'Description'].decode())
            return match.group(1).split('|')

def parse_csq_entries(csq_entries, csq_fields):
    transcripts = {}
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_fields, values):
            transcript[key] = value
        if transcript['Allele'] not in transcripts.keys():
            transcripts[transcript['Allele']] = []
        transcripts[transcript['Allele']].append(transcript)
    return transcripts

def resolve_alleles(entry, csq_alleles):
    alleles = {}
    if entry.is_indel:
        for alt in entry.ALT:
            alt = str(alt)
            if alt[0:1] != entry.REF[0:1]:
                csq_allele = alt
            elif alt[1:] == "":
                csq_allele = '-'
            else:
                csq_allele = alt[1:]
            alleles[alt] = csq_allele
    elif entry.is_sv:
        for alt in alts:
            if len(alt) > len(entry.REF) and 'insertion' in csq_alleles:
                alleles[alt] = 'insertion'
            elif len(alt) < len(entry.REF) and 'deletion' in csq_alleles:
                alleles[alt] = 'deletion'
            elif len(csq_alleles) == 1:
                alleles[alt] = list(csq_alleles)[0]
    else:
        for alt in entry.ALT:
            alt = str(alt)
            alleles[alt] = alt
    return alleles

def transcript_for_alt(transcripts, alt):
    for transcript in transcripts[alt]:
        if transcript['PICK'] == '1':
            return transcript
    return transcripts[alt][0]

def decode_hex(string):
    hex_string = string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

(script, tsv_filename, vcf_filename, vep_fields, output_dir) = sys.argv
vep_fields_list = vep_fields.split(',')

vcf_file = VCF(vcf_filename)

csq_fields = parse_csq_header(vcf_file)

vep = {}
for variant in vcf_file:
    chr = str(variant.CHROM)
    pos = str(variant.POS)
    ref = str(variant.REF)
    alts = variant.ALT

    if chr not in vep:
        vep[chr] = {}

    if pos not in vep[chr]:
        vep[chr][pos] = {}

    if ref not in vep[chr][pos]:
        vep[chr][pos][ref] = {}

    csq = variant.INFO.get('CSQ')
    if csq is not None:
        transcripts = parse_csq_entries(csq.split(','), csq_fields)
    alleles_dict = resolve_alleles(variant, transcripts.keys())
    for alt in alts:
        if alt not in vep[chr][pos][ref]:
            if transcripts is not None and alleles_dict[alt] in transcripts:
                vep[chr][pos][ref][alt] = transcript_for_alt(transcripts, alleles_dict[alt])
            else:
                vep[chr][pos][ref][alt] = None
        else:
            sys.exit("VEP entry for at CHR %s, POS %s, REF %s , ALT % already exists" % (chr, pos, ref, alt) )


with open(tsv_filename, 'r') as input_filehandle:
    reader = csv.DictReader(input_filehandle, delimiter = "\t")
    output_filehandle = open(os.path.join(output_dir, 'variants.annotated.tsv'), 'w')
    writer = csv.DictWriter(output_filehandle, fieldnames = reader.fieldnames + vep_fields_list, delimiter = "\t")
    writer.writeheader()
    for entry in reader:
        row = entry
        for field in vep_fields_list:
            field_annotations = []
            for alt in entry['ALT'].split(','):
                vep_annotations = vep[entry['CHROM']][entry['POS']][entry['REF']][alt]
                if vep_annotations is not None and field in vep_annotations:
                    annotation = vep_annotations[field]
                    decoded_annotation = re.sub(r'%[0-9|A-F][0-9|A-F]', decode_hex, annotation)
                    field_annotations.append(decoded_annotation)
                else:
                    field_annotations.append('-')
            row[field] = ','.join(field_annotations)
        writer.writerow(row)
    output_filehandle.close()
