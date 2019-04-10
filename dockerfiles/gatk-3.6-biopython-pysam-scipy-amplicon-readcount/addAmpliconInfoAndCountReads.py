from __future__ import division
import pysam
from Bio import pairwise2
import argparse
from string import maketrans
from scipy.stats import binom
import re

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('bamfile',help='Consensus BAM file')
parser.add_argument('samplename',help='Sample name')
parser.add_argument('-r',"--reference",default="/gscuser/dspencer/refdata/GRCh37/all_sequences.fa",
                    help='reference genome fasta file')
parser.add_argument('-w',"--window",type=int,default=150,
                    help='window for creating ref and alt sequences')
parser.add_argument('-m',"--minreads",type=int,default=5,
                    help='minimum reads for defining an amplicon')
parser.add_argument('-v',"--minreads4vaf",type=int,default=5,
                    help='minimum reads in an amplicon for amplicon-based vaf calculation')
parser.add_argument('-a',"--minampnumber",type=int,default=1,
                                        help='Minimum number of covering amplicons with support for variant allele')
parser.add_argument('-n',"--maxreadnm",type=int,default=4,
                                        help='Maximum edit distance from NM tag for a read to be counted')
parser.add_argument('-q',"--minbasequal",type=int,default=15,
                                        help='Minium base quality to assess variant')
parser.add_argument('-p',"--strandpvalue",type=float,default=0.05,
                                        help='Binomial p-value for strand bias')
parser.add_argument('-s',"--minstrandreads",type=int,default=5,
                                        help='Mininum reads on each strand if sufficient strand coverage is present')
parser.add_argument('-d',"--minvaf",type=float,default=0.02,
                                        help='Mininum VAF required')

args = parser.parse_args()

mysample = args.samplename

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)
# open bam file
samfile = pysam.AlignmentFile(args.bamfile,"rb")
# open reffasta
fa = pysam.FastaFile(args.reference)

# window on either side of a variant for indel annotation 
window = args.window
# minimum number of read pairs with a unique insert length to define an amplicon
minampreads = args.minreads
# minimum number of reads in an amplicon needed to be eligible to calculate a VAF
minampreadsforvaf = args.minreads4vaf

minampsupport = args.minampnumber

maxNM = args.maxreadnm

minqual = args.minbasequal

minmapqual = 1

strandp = args.strandpvalue

minstrandsupport = args.minstrandreads

minvaf = args.minvaf

vcffile.header.filters.add("AMPSupport",None,None,'Fails requirement of having >='+str(minampsupport)+' amplicons with support for the variant')
vcffile.header.filters.add("LowReads",None,None,'Less than '+str(minampreads)+' Q'+str(minqual)+' bases with variant allele among reads with NM < '+str(maxNM))
vcffile.header.filters.add("StrandSupport",None,None,'Variant allele support is lacking on one strand despite adequate coverage (binomial P < '+str(strandp)+' for observing at least '+str(minampreads)+' given the coverage on that strand')
vcffile.header.formats.add("TAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position')
vcffile.header.formats.add("SAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position that support the alternate allele')
vcffile.header.formats.add("AMPS", 1, 'String', 'Amplicon string, in the format amplicon id 1,num. reference alleles in amplicon 1, num. alternate alleles in amplicon 1;amplicon id 2,num ref, num alt;...')
vcffile.header.formats.add("CVAF", 1, 'Float', 'Calculated VAF, which is based only on supporting amplicons with with >' + str(minampreadsforvaf) + ' reads')
vcffile.header.add_line("##ampliconparseroptions={'minampfraction':"+str(minampsupport)+"'window':"+str(window)+",'minampreads':"+str(minampreads)+",'minampreadsforvaf':"+str(minampreadsforvaf)+"}")

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print "\n".join(hdr) + "\n" + ("\t".join((hd.split("\t"))[0:9])) + "\t" + mysample

for vline in vcffile.fetch():

    for alt in vline.alts:

        rec = vline
        rec.alts = [ alt ]
    
        amps = {'ref':[],'alt':[]}
        ampcnt = {}
        amptcnt = {}
        ampstring = list()
        ampnum = 0
        ampvafs = list()
        strands = {'ref':[],'alt':[]}
        
        # if the variant is a substitution
        if len(rec.ref) == len(rec.alts[0]) and len(rec.ref) == 1:
            
            for pileup in samfile.pileup(contig=rec.contig, start=rec.pos-1, stop=rec.pos):
                if pileup.pos == rec.pos-1: # only consider the variant position
                    
                    for read in pileup.pileups:

                        # skip if more than maxNM edit distance for this read or position in indel or 
                        if read.alignment.get_tag("NM") > maxNM + 1 or read.is_del or read.is_refskip or int(read.alignment.query_qualities[read.query_position]) < minqual or read.alignment.mapping_quality < minmapqual:
                            continue
                        
                        # get number of raw reads for this collapsed read family
                        numbcreads = read.alignment.get_tag("BC").split(":")[0]
                        
                        # get the amplicon using the template length and read 1 strand
                        ampname = read.alignment.get_tag("X1")
                        
                        # count reads per amplicon
                        if ampname not in ampcnt:
                            ampcnt[ampname] = 1
                        else:
                            ampcnt[ampname] += 1
                            
                        # count reads per amplicon
                        if ampname not in amptcnt:
                            amptcnt[ampname] = numbcreads
                        else:
                            amptcnt[ampname] += numbcreads
                            
                        # count alleles by amplicon
                        if read.alignment.query_sequence[read.query_position] != rec.alts[0]:
                            amps['ref'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['ref'].append(1)
                            else:
                                strands['ref'].append(0)
                                
                        elif read.alignment.query_sequence[read.query_position] == rec.alts[0]:
                            amps['alt'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['alt'].append(1)
                            else:
                                strands['alt'].append(0)
                            
                        else:
                            amps['ref'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['ref'].append(1)
                            else:
                                strands['ref'].append(0)
                                                
        else: # if indel
                        
            refseqstart = rec.pos-window-1
            refseqend = rec.pos+window
            refseq = fa.fetch(rec.contig,refseqstart,refseqend)
            altseq = fa.fetch(rec.contig,refseqstart,rec.pos-1) # 0-based
                
            varlen = 0
                
            # if deletion
            if len(rec.ref) > len(rec.alts[0]):
                altseq = altseq + fa.fetch(rec.contig,rec.pos+len(rec.ref)-1,rec.pos+window+len(rec.ref))
                varlen = len(rec.ref) - len(rec.alts[0])
                    
            else: # if insertion
                altseq = altseq + rec.alts[0] + fa.fetch(rec.contig,rec.pos,rec.pos+window)
                varlen = len(rec.alts[0]) - len(rec.ref)
                    
            reads = {}
            for pileup in samfile.pileup(contig=rec.contig, start=rec.pos-1, stop=rec.pos):
                if pileup.pos == rec.pos-1:
                    for read in pileup.pileups:

                        # skip low map quality reads
                        if read.alignment.mapping_quality < minmapqual:
                            continue
                        
                        # skip if more than maxNM + variant length edit distance for this read                                                                                                                                               
                        if read.alignment.get_tag("NM") > maxNM + varlen:
                            continue
                                        
                        # get number of raw reads for this collapsed read family
                        numbcreads = read.alignment.get_tag("BC").split(":")[0]
                        
                        # get the amplicon using the template length and read 1 strand
                        ampname = read.alignment.get_tag("X1")
                                
                        # count reads per amplicon
                        if ampname not in ampcnt:
                            ampcnt[ampname] = 1
                        else:
                            ampcnt[ampname] += 1
                            
                        # count reads per amplicon
                        if ampname not in amptcnt:
                            amptcnt[ampname] = numbcreads
                        else:
                            amptcnt[ampname] += numbcreads
                                                
                        ampname = read.alignment.get_tag("X1")
                                
                        # if the read has fewer mismatches (NM tag) than the indel length and has no softclipped bases
                        # then assign to reference allele
                        if len(read.alignment.cigartuples) == 1 and read.alignment.cigartuples[0][0] == 0 and read.alignment.cigartuples[0][1] == read.alignment.query_length and read.alignment.get_cigar_stats()[0][10] < abs(len(rec.ref) - len(rec.alts[0])):
                            amps['ref'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['ref'].append(1)
                            else:
                                strands['ref'].append(0)
                                        
                            # if the record is a deletion and the read has a deletion at this position
                            # with the proper length then assign to the alt allele
                        elif len(rec.ref) > len(rec.alts[0]) and read.indel < 0 and abs(read.indel) == len(rec.ref)-len(rec.alts[0]):
                            amps['alt'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['alt'].append(1)
                            else:
                                strands['alt'].append(0)                        
                                        
                        # if the record is an insertion and the read has an insertion
                        # with the proper bases and length then assign to the alt allele
                        elif len(rec.ref) < len(rec.alts[0]) and read.indel > 0 and rec.alts[0] == read.alignment.seq[read.query_position:read.query_position+read.indel+1]:
                            amps['alt'].append(ampname)
                            if read.alignment.is_reverse:
                                strands['alt'].append(1)
                            else:
                                strands['alt'].append(0)
                                        
                        # if none of the above are satisified, use the alignment method
                        # to assign reads to alleles
                        else:
                            rdseq = read.alignment.query_sequence
                            
                            alnref = pairwise2.align.localms(rdseq,refseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                            alnalt = pairwise2.align.localms(rdseq,altseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                                    
                            if alnref >= alnalt:
                                amps['ref'].append(ampname)
                                if read.alignment.is_reverse:
                                    strands['ref'].append(1)
                                else:
                                    strands['ref'].append(0)
                                            
                            else:
                                amps['alt'].append(ampname)
                                if read.alignment.is_reverse:
                                    strands['alt'].append(1)
                                else:
                                    strands['alt'].append(0)
                                            
        # count number of reads per amplicon and stor per amplicon data if there are enough
            
        ampwithalt = 0
        nr=0
        nv=0
        fwdstrandvaf=0.0
        revstrandvaf=0.0
        myvaf = 0.0
        
        ampnum = len(ampcnt)
            
#        if(ampnum < 1):
#            sys.exit('No amplicons detected')
            
        for a in ampcnt:
            if ampcnt[a] >= minampreads:
                if amps['alt'].count(a) > 0:
                    ampwithalt += 1
                    nr += amps['ref'].count(a)
                    nv += amps['alt'].count(a)                            
                    ampstring.append(str(a) + "," + str(amps['ref'].count(a)) + "," + str(amps['alt'].count(a)))
                                                                            
        if (nr+nv) > 0:
            myvaf = nv / (nv + nr)
                
        if strands['alt'].count(0) + strands['ref'].count(0) > 0:
            fwdstrandvaf = strands['alt'].count(0) / (strands['alt'].count(0) + strands['ref'].count(0))
                                                                                    
        if strands['alt'].count(1) + strands['ref'].count(1) > 0:        
            revstrandvaf = strands['alt'].count(1) / (strands['alt'].count(1) + strands['ref'].count(1))
                                                                                        
        mygt = ('.','.')
    
        # filter variant, if necessary
                                                                                        
        # special case for PindelITD
        isitd = str(rec.info["set"]).find("ITD")
    
        # if its called by PindelITD then only require at least 1 read on each strand
        if isitd > 0 and len(rec.alts[0]) > len(rec.ref) and strands['alt'].count(0) > 0 and strands['alt'].count(1) > 0:
            rec.filter.clear()
            rec.filter.add("PASS")
        
        # minimum number of supporting amplicons
        elif ampnum < minampsupport:
            rec.filter.add("AMPSupport")
        
        # min vaf or variant supporting read count
        elif myvaf < minvaf or nv < minampreads:
            rec.filter.add("LowReads")

        # minimum strand support: at least 5 reads on each strand, unless there are too few reads on that strand to expect 5 supporting read (via binomial prob.) 
        elif ((strands['alt'].count(1) < minstrandsupport and binom.cdf(minstrandsupport, strands['alt'].count(1)+strands['ref'].count(1),fwdstrandvaf, loc=0) < strandp) or (strands['alt'].count(0) < minstrandsupport and binom.cdf(minstrandsupport, strands['alt'].count(0)+strands['ref'].count(0),revstrandvaf, loc=0) < strandp)):
            rec.filter.add("StrandSupport")

        else:
            rec.filter.clear()
            rec.filter.add("PASS")
            
        for s in rec.samples:
            if rec.samples[s]['GT'] == (1,1) and myvaf > .99:
                mygt = (1,1)
            else:
                mygt = (0,1)
            
        mysample=0
        rec.samples[mysample]['ADF'] = strands['alt'].count(0)
        rec.samples[mysample]['ADR'] = strands['alt'].count(1)
        rec.samples[mysample]['RDF'] = strands['ref'].count(0)
        rec.samples[mysample]['RDR'] = strands['ref'].count(1)
        rec.samples[mysample]['GT'] = mygt
        rec.samples[mysample]['NR'] = len(amps['alt']) + len(amps['ref'])
        rec.samples[mysample]['NV'] = len(amps['alt'])
        rec.samples[mysample]['TAMP'] = ampnum
        rec.samples[mysample]['SAMP'] = ampwithalt
        rec.samples[mysample]['AMPS'] = ';'.join(ampstring)
        rec.samples[mysample]['CVAF'] = myvaf

        print "\t".join(str(rec).rstrip().split("\t")[0:10])

    # end of for each alts
    
# end vcf fetch

fa.close()
vcffile.close()
samfile.close()
