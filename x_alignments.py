##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

import os
import sys
from subprocess import *
from multiprocessing import Pool
import pysam


#
def unwrap_self_extract_reads_for_region(arg, **kwarg):
    return BamInfo.extract_reads_for_region(*arg, **kwarg)

def unwrap_self_extract_mate_reads_for_region(arg, **kwarg):
    return BamInfo.extract_mate_reads_of_region(*arg, **kwarg)

####Function: Functions for processing alignment, like trim, change format, and etc.
class BamInfo():#
    def __init__(self, sf_bam, sf_ref):
        self.sf_bam = sf_bam
        self.out_header = None
        self.chrm_id_name = {}
        self.sf_reference=sf_ref
        # if the chromosome name in format: chr1, then return true,

    def index_reference_name_id(self):
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        header = bamfile.header
        bamfile.close()

        l_chrms = header['SQ']
        chrm_id = 0
        for record in l_chrms:
            chrm_name = record['SN']
            self.chrm_id_name[chrm_id] = chrm_name
            chrm_id += 1

    def get_all_reference_names(self):
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        header = bamfile.header
        bamfile.close()

        l_chrms = header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_name = record['SN']
            m_chrms[chrm_name] = 1
        return m_chrms
####
    # get the chrom name and length in a dictionary
    def get_all_chrom_name_length(self):
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        header = bamfile.header
        bamfile.close()

        l_chrms = header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_name = record['SN']
            chrm_length = int(record['LN'])
            m_chrms[chrm_name] = chrm_length
        return m_chrms

    def get_reference_name_by_id(self, id):#
        i_id = int(id)
        if i_id in self.chrm_id_name:
            return self.chrm_id_name[i_id]
        else:
            return "*"

####
    # else if in format: 1, then return false
    def is_chrm_contain_chr(self):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        self.out_header = samfile.header  # it's a dictionary
        samfile.close()

        l_chrms = self.out_header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_id = record['SN']
            m_chrms[chrm_id] = 1

        b_with_chr = False
        if "chr1" in m_chrms:
            b_with_chr = True

        self.b_with_chr = b_with_chr
        return b_with_chr




    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr #######################################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    # calculate coverage
    def calc_coverage(self, fpath, read_length, genome_length):
        cnt_lines = 0
        cmd = ""
        if fpath.lower().endswith(('.fq', '.fastq')):
            cmd = "echo $(wc -l {0})".format(fpath)
        elif fpath.lower().endswith(('.fastq.gz', 'fq.gz', 'gz')):
            cmd = "echo $(zcat {0} | wc -l)".format(fpath)
        else:
            print(("Something wrong with the raw reads files format:", fpath))

        tp = tuple(Popen(cmd, shell=True, stdout=PIPE).communicate())
        lcnt = str(tp[0]).split()

        cnt_lines = int(lcnt[0])
        cnt_reads = int(cnt_lines) / 4
        cov = float(cnt_reads * read_length) / float(genome_length)
        return cov

    ##read in the read names from file
    def load_read_names(self, sf_names):
        with open(sf_names, 'r') as infile:
            l_names = infile.read().splitlines()
        if '' in l_names:
            l_names.remove('')
        return l_names

    # extract the reads by read names
    def extract_reads_by_name_list(self, sf_names, sf_bam, sf_out_bam):
        l_names = self.load_read_names(sf_names)
        bamfile = pysam.AlignmentFile(sf_bam, 'rb', reference_filename=self.sf_reference)
        name_indexed = pysam.IndexedReads(bamfile)  # here use hashing to save the read names in the memory
        name_indexed.build()
        header = bamfile.header.copy()
        out = pysam.Samfile(sf_out_bam, 'wb', header=header)
        for name in l_names:
            try:
                name_indexed.find(name)
            except KeyError:
                pass
            else:
                iterator = name_indexed.find(name)
                for x in iterator:  # x is an alignment
                    ###here need to check whether this is the first or second read we wanted!!!!
                    out.write(x)
        out.close()

    def extract_reads_for_region(self, record):
        chrm=record[0]
        istart=record[1]
        iend=record[2]
        sf_read_names=record[3]
        s_working_folder=record[4]
        if s_working_folder[-1]!="/":
            s_working_folder+="/"

        #first write the bam to file, and then index it, and finally get the fastq reads
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        sf_tmp_bam=s_working_folder+"{0}_{1}_{2}.bam".format(chrm, istart, iend)
        out_tmp_bam = pysam.Samfile(sf_tmp_bam, 'wb', template=bamfile)
        for algnmt in bamfile.fetch(chrm, istart, iend):
            out_tmp_bam.write(algnmt)
        out_tmp_bam.close()
        bamfile.close()

        sf_parsed_bam=s_working_folder+"{0}_{1}_{2}.disc.bam".format(chrm, istart, iend)
        self.extract_reads_by_name_list(sf_read_names, sf_tmp_bam, sf_parsed_bam)

    #convert an alignmet to a fasta format seq
    def cvt_alignmt_to_fa(self, alignmt):
        read_name=alignmt.query_name
        query_seq = alignmt.query_sequence

####

    ####give the read name file, return the read fall in the given region
    def _parse_read_name_pos_of_region(self, sf_name_pos, chrm, start, end):
        m_reads={}
        with open(sf_name_pos) as fin_name_pos:
            #each line in format: chrm, map_pos, is_rc, is_mate_rc, mate_chrm, mate_pos, query_name, s_mate_rc, insertion_pos
            for line in fin_name_pos:
                fields=line.split()
                ins_chrm = fields[0]
                anchor_read_map_pos=int(fields[1]) ##this is the map_pos of the anchor read
                s_anchor_rc=fields[2] # "0" or "1"
                s_anchor_mate_rc=fields[3]  # "0" or "1"
                mate_chrm=fields[4]
                mate_pos=int(fields[5])

                if chrm != mate_chrm:
                    continue
                if mate_pos<(start-100) or mate_pos>(end+100):
                    continue

                rname=fields[6]
                s_first=fields[7]
                s_insertion_pos=fields[8]
                b_lclip=fields[9] # 0 or 1
                b_rclip=fields[10] # 0 or 1

                s_insertion_site = "{0}~{1}".format(ins_chrm, s_insertion_pos)
                s_info="{0}~{1}".format(rname, s_first)
                if s_info not in m_reads:
                    m_reads[s_info]={}
                m_reads[s_info][s_insertion_site]=(anchor_read_map_pos, s_anchor_rc, s_anchor_mate_rc, b_lclip, b_rclip)
        return m_reads
#

    def extract_mate_reads_of_region(self, record):
        chrm=record[0]
        bin_start=int(record[1])
        bin_end=int(record[2])
        sf_name_pos=record[3]
        s_working_folder=record[4]

        m_reads=self._parse_read_name_pos_of_region(sf_name_pos, chrm, bin_start, bin_end)

        #first write the bam to file, and then index it, and finally get the fastq reads
        bamfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        sf_tmp_disc_fa=s_working_folder+"{0}_{1}_{2}.tmp.disc.fa".format(chrm, bin_start, bin_end)
        with open(sf_tmp_disc_fa,"w") as fout_fa:
            for alignmt in bamfile.fetch(chrm, bin_start, bin_end):
                s_read_name = alignmt.query_name
                s_map_pos= str(alignmt.reference_start)
                read_seq = alignmt.query_sequence
                is_first = 0
                if alignmt.is_read1 == True:
                    is_first = 1

                s_read_id = "{0}~{1}".format(s_read_name, is_first)
                if s_read_id not in m_reads:
                    continue
                for s_insertion in m_reads[s_read_id]:
                    (anchor_map_pos, s_anchor_rc, s_anchor_mate_rc, b_lclip, b_rclip)=m_reads[s_read_id][s_insertion]
                    s_read_info = ">{0}~{1}~{2}~{3}~{4}~{5}~{6}\n".format(s_read_id, b_lclip, b_rclip, s_anchor_rc,
                                                                          s_anchor_mate_rc, anchor_map_pos, s_insertion)
                    fout_fa.write(s_read_info)
                    fout_fa.write(read_seq + "\n")
        bamfile.close()
####

