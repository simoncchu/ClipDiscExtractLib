import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_alignments import *
from intervaltree import *
import global_values


def unwrap_self_collect_clip_disc_reads(arg, **kwarg):
    return TClipDisc.collect_clipped_disc_reads_by_region(*arg, **kwarg)

####
class TClipDisc():####
    def __init__(self, sf_bam, working_folder, n_jobs, sf_ref):
        self.sf_bam = sf_bam
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.sf_reference = sf_ref
        self.l_allowed_rep=["SVA","HERV","LTR", "LINE1", "L1"]#in consensus algnmnt, decoy sequence are added, if algined to them, then filter out
        self.MIN_MAPQ_RNA=50

    def collect_clipped_disc_reads_of_given_list(self, l_chrm_records, bin_size, sf_all_clip_fq, sf_disc_fa):
        pool = Pool(self.n_jobs)
        l_low_mapq = pool.map(unwrap_self_collect_clip_disc_reads,
                              list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()
        #for rcd in l_chrm_records:
        #    self.collect_clipped_disc_reads_by_region(rcd)
####
        ###II. For disc reads
        ### here merge all the read_id files and pos files
        ####in format:chrm, map_pos, mate_chrm, mate_pos, query_name, s_mate_first, insertion_pos
        ####it is possible there are duplicate records are merged
        sf_all_disc_pos = self.working_folder + "all_disc_pos" + global_values.DISC_POS_SUFFIX  # this is to save all the disc pos
        with open(sf_all_disc_pos, "w") as fout_merged_disc_pos:
            # m_read_names={}####in case duplicate reads are saved
            for record in l_chrm_records:
                chrm = record[0]
                insertion_pos = record[1]
                working_folder= record[8]
                s_pos_info = "{0}_{1}".format(chrm, insertion_pos)
                sf_disc_pos = working_folder + s_pos_info + global_values.DISC_POS_SUFFIX  # this is to save the disc positions
                if os.path.isfile(sf_disc_pos) == True:
                    with open(sf_disc_pos) as fin_disc_pos:
                        for line in fin_disc_pos:
                            fout_merged_disc_pos.write(line)
                if os.path.isfile(sf_disc_pos):
                    os.remove(sf_disc_pos)

        # now, need to retrieve the reads according to the read names and disc positions
        bam_info = BamInfo(self.sf_bam, self.sf_reference)
        bam_info.extract_mate_reads_by_name(sf_all_disc_pos, bin_size, self.working_folder, self.n_jobs, sf_disc_fa)
        return l_low_mapq
####

    ###Problem here: 1. chrm in "candidate_list" may not consistent with chrm in bam file
    ###2. all should follow the style in candidate list
    def collect_clipped_disc_reads_by_region(self, record):##
        #record in format: (s_chrm, i_event_start, i_event_end, i_left_most, i_right_most, l_left_rgns, l_right_rgns,
                            #sf_algnmt, sf_new_wfolder, sf_black_list)
        chrm = record[0]  ##this is the chrm style in candidate list
        insertion_pos = record[1] #or this is the left breakpoint for a reference copy
        i_copy_end=record[2]
        b_polymorphic=False
        if i_copy_end==-1:
            b_polymorphic=True
        start_pos = record[3]
        end_pos = record[4] #
        l_left_rgns=record[5]
        l_right_rgns=record[6]
        sf_bam = record[7]
        working_folder = record[8] #
        sf_black_list = record[9]#

        interval_tree = IntervalTree()
        for (i_tmp_start, i_tmp_end) in l_left_rgns:
            interval_tree.addi(i_tmp_start-1, i_tmp_end+1)
        for (i_tmp_start, i_tmp_end) in l_right_rgns:
            interval_tree.addi(i_tmp_start-1, i_tmp_end+1)

        bam_info = BamInfo(sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        chrm_in_bam = bam_info.process_chrm_name(chrm, b_with_chr)

        # load the reads, and write the related clipped part into file
        s_pos_info = "{0}_{1}".format(chrm, insertion_pos)
        sf_disc_pos = working_folder + s_pos_info + global_values.DISC_POS_SUFFIX  # this is to save the discordant positions
        f_disc_pos = open(sf_disc_pos, "w")


        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        m_chrm_id=self._get_chrm_id_name(samfile)##
        i_max_rlth=0
        for algnmt in samfile.fetch(chrm_in_bam, start_pos, end_pos):  ##fetch reads mapped to "chrm:start_pos-end_pos"
            ##here need to skip the secondary and supplementary alignments?
            # if algnmt.is_secondary or algnmt.is_supplementary:
            #     continue
            #if algnmt.is_duplicate == True:  ##duplciate
            #    continue
            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            if algnmt.is_unmapped == True:  # unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue

            ####
            if algnmt.mapping_quality < self.MIN_MAPQ_RNA:
                continue

            s_read_seq = algnmt.query_sequence
            i_rlen=len(s_read_seq)
            if i_rlen > i_max_rlth:
                i_max_rlth = i_rlen

            mate_chrm = '*'
            mate_pos = 0
            if (algnmt.next_reference_id in m_chrm_id) and (algnmt.mate_is_unmapped == False) \
                    and (algnmt.next_reference_id >= 0):
                mate_chrm = algnmt.next_reference_name
                mate_pos = algnmt.next_reference_start
            if mate_chrm == "*":  ##unmapped reads are not interested!
                continue

            query_name = algnmt.query_name
            #Here we require at least half of the read is overlapped with the exon
            map_pos = algnmt.reference_start
            b_fully_mapped = False
            i_half_read=int(i_max_rlth/2)
            i_check_start = map_pos + i_half_read
            l_hits=interval_tree[i_check_start] #if outside the selected regions, then skip
            if len(l_hits)<=0:#head of read is not in range
                i_map_end = -1
                if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                    b_fully_mapped = True
                    i_map_end = map_pos + l_cigar[0][1]
                else:
                    for t_cigar_rcd in l_cigar:
                        i_flag = t_cigar_rcd[0]
                        if i_flag != 2:
                            i_map_end += t_cigar_rcd[1]
                i_check_end=i_map_end-i_half_read
                l_end_hits=interval_tree[i_check_end]
                if len(l_end_hits)<=0:#end also not in range
                    continue
####
            is_rc = 0
            if algnmt.is_reverse == True:  # is reverse complementary
                is_rc = 1
####
            b_left_clip = 0  # whether is left clip
            if l_cigar[0][0] == 4:  # left clipped
                b_left_clip = 1
            b_right_clip=0#whether is right clip
            if l_cigar[-1][0] == 4:  # right clipped
                b_right_clip = 1

            # we need check "is_mate_rc", because reads from alignment is always same as reference
            # so need this information to get the original orientation
            is_mate_rc = 0
            if algnmt.mate_is_reverse == True:  #mate is reverse complementary
                is_mate_rc = 1
            ## here only collect the read names for discordant reads, later will re-align the discordant reads
            i_dist_allow_any=10
            if self.is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos, i_dist_allow_any) == True:
                # check where the mate is mapped, if within a repeat copy, then get the position on consensus
                # f_disc_names.write(query_name + "\n")
                s_mate_first = 1  # whether the mate read is the "first read" in a pair
                if b_first == True:
                    s_mate_first = 0

                # here mate_chrm must be the style in the bam file
                # And chrm must be the style in the candidate file
                s_mate_pos_info = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}" \
                                  "\t{9}\t{10}\n".format(chrm, map_pos, is_rc, is_mate_rc, mate_chrm,
                                                         mate_pos, query_name, s_mate_first, insertion_pos,
                                                         b_left_clip, b_right_clip)
                f_disc_pos.write(s_mate_pos_info)
        samfile.close()
        f_disc_pos.close()
        return (chrm, insertion_pos)
####
####
    # For this version, reads are aligned to the consensus sequence
    # Note that the poly-A part in the consensus has been removed, so it's possible reads will be clipped mapped (at end)
    def parse_disc_algnmt_consensus(self, sf_disc_alignmt, bmapped_cutoff):#
        samfile = pysam.AlignmentFile(sf_disc_alignmt, "r", reference_filename=self.sf_reference)
        m_disc_pos = {}
        for algnmt in samfile.fetch():  #
            # mapq = algnmt.mapping_quality
            # if mapq<global_values.MINIMUM_DISC_MAPQ:##############Here should be very careful for SVA and Alu!!!!!!!!!!!!!!!!!!!!!!!
            #     continue
            ####also, for clipped mapped reads, need to check the clipped parts whether can be split to two parts!!!!!!
            # fmt:read_id~is_first~is_anchor_rc~is_anchor_mate_rc~anchor_pos~s_insertion_chrm~s_insertion_pos
            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue

            b_allowed = False
            s_hit_rep_name=algnmt.reference_name
            for s_allowed_rep in self.l_allowed_rep:
                if s_allowed_rep in s_hit_rep_name:
                    b_allowed=True
                    break
            if b_allowed==False:
                continue

            # first check whether read is qualified mapped
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)
            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            read_info = algnmt.query_name
            read_info_fields = read_info.split(global_values.SEPERATOR)
            s_anchor_lclip = read_info_fields[-7]  # anchor read is left clip or not: 1 indicates clip
            s_anchor_rclip = read_info_fields[-6]  # anchor read is right clip or not: 1 indicates clip
            s_anchor_rc = read_info_fields[-5]
            s_anchor_mate_rc = read_info_fields[-4]
            anchor_map_pos = int(read_info_fields[-3])  # this is the position of the left-most mapped base
            ins_chrm = read_info_fields[-2]
            ins_pos = int(read_info_fields[-1])
            read_seq = algnmt.query_sequence

            # map position on repeat consensus
            pos_in_consensus = algnmt.reference_start + len(read_seq) / 2
            if ins_chrm not in m_disc_pos:
                m_disc_pos[ins_chrm] = {}
            if ins_pos not in m_disc_pos[ins_chrm]:
                m_disc_pos[ins_chrm][ins_pos] = []

            ##Here to check orientation of the two reads, and later will be used to detect inversion
            b_anchor_mate_rc = False
            if s_anchor_mate_rc == "1":
                b_anchor_mate_rc = True
            b_masked_rc = algnmt.is_reverse
            if b_anchor_mate_rc == True:
                b_masked_rc = (not algnmt.is_reverse)

            b_anchor_rc = False
            if s_anchor_rc == "1":  ####
                b_anchor_rc = True

            b_anchor_lclip = False
            if s_anchor_lclip == "1":
                b_anchor_lclip = True

            b_anchor_rclip = False
            if s_anchor_rclip == "1":
                b_anchor_rclip = True

            b_same_dir_xor = (b_masked_rc != b_anchor_rc)  # check the xor results
            m_disc_pos[ins_chrm][ins_pos].append((anchor_map_pos, b_anchor_lclip, b_anchor_rclip, pos_in_consensus,
                                                  b_same_dir_xor, b_anchor_rc))

        samfile.close()
        return m_disc_pos
####
####
####
####
    def is_two_side_clipped(self, l_cigar, i_min_clip):
        b_two_clip=False
        if l_cigar[-1][0] == 4 and l_cigar[0][0] == 4:#both side clip
            if l_cigar[-1][1] > i_min_clip and l_cigar[0][1] > i_min_clip:
                return True
        return b_two_clip
    ####
    ##whether a pair of read is discordant (for TEI only) or not
    def is_discordant(self, chrm, map_pos, mate_chrm, mate_pos, is_threshold):
        # b_disc=False
        if chrm != mate_chrm:  ###of different chroms
            return True
        else:
            if abs(mate_pos - map_pos) > is_threshold:  # Of same chrom, but insert size are quite large
                return True
                # if first_rc==second_rc: ##direction are abnormal, by default consider (F,R) as right mode
                #     return True
        return False
####
####
    def _cvt_to_Ascii_quality(self, l_score):
        new_score = [x + 33 for x in l_score]
        return ''.join(map(chr, new_score))

    def _is_qualified_clip(self, l_score):
        n_char=len(l_score)
        n_half=n_char/2
        n_cnt=0
        for i_score in l_score:
            if i_score < global_values.CLIP_PHRED_SCORE_CUTOFF:
                n_cnt+=1
            if n_cnt>n_half:
                return False
        return True
    def _get_chrm_id_name(self, samfile):
        m_chrm = {}
        references = samfile.references
        for schrm in references:
            chrm_id = samfile.get_tid(schrm)
            m_chrm[chrm_id] = schrm
        m_chrm[-1] = "*"
        return m_chrm

    ####
    # check the clipped part is qualified aligned or not
    def is_clipped_part_qualified_algnmt(self, l_cigar, ratio_cutoff):
        if len(l_cigar) < 1:  # wrong alignment
            return False, 0
        if len(l_cigar) > 2:
            ####check the cigar
            ###if both clipped, and the clipped part is large, then skip
            b_left_clip = False
            i_left_clip_len = 0
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                b_left_clip = True
                i_left_clip_len = l_cigar[0][1]
            b_right_clip = False
            i_right_clip_len = 0
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                b_right_clip = True
                i_right_clip_len = l_cigar[-1][1]

            if b_left_clip == True and b_right_clip == True:
                if (i_left_clip_len > global_values.MAX_CLIP_CLIP_LEN) and (
                        i_right_clip_len > global_values.MAX_CLIP_CLIP_LEN):
                    return False, 0

        ####for the alignment (of the clipped read), if the mapped part is smaller than the clipped part,
        ####then skip
        n_total = 0
        n_map = 0
        for (type, lenth) in l_cigar:
            if type == 0:
                n_map += lenth
            if type != 2:  # deletion is not added to the total length
                n_total += lenth

        if n_map < (n_total * ratio_cutoff):###########require at least 3/4 of the seq is mapped !!!!!!!!
            return False, 0
        return True, n_map
####