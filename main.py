# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import sys
from clip_disc_collection import *

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


def _prepare_records(self, l_focal_regions, sf_algnmt, sf_black_list):
    l_records=[]
    sf_new_wfolder=self.s_wfolder +"tmp/"
    if os.path.exists(sf_new_wfolder)==False:
        os.mkdir(sf_new_wfolder)
    for tmp_rcd in l_focal_regions:
        s_chrm=tmp_rcd[0]
        i_event_start=tmp_rcd[1]
        i_event_end=tmp_rcd[2]
        l_left_rgns=tmp_rcd[3]
        l_right_rgns = tmp_rcd[4]
        #print(l_left_rgns)

        i_left_most=self._get_left_most_position(l_left_rgns)
        i_right_most=self._get_right_most_position(l_right_rgns)
        new_rcd=(s_chrm, i_event_start, i_event_end, i_left_most, i_right_most, l_left_rgns, l_right_rgns,
                 sf_algnmt, sf_new_wfolder, sf_black_list)
        l_records.append(new_rcd)
    return l_records

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

    sf_rna_algnmt=sys.argv[1]
    s_wfolder=sys.argv[2]
    n_jobs=int(sys.argv[3])
    sf_ref=sys.argv[4]

    sf_black_list=None
    l_all_slct_regions=[]
    l_left_slcts=[]
    l_right_slcts=[]
    l_all_slct_regions.append((s_chrm, i_event_start, i_event_end, l_left_slcts, l_right_slcts))
    l_record=_prepare_records(l_all_slct_regions,sf_rna_algnmt, sf_black_list)


    xclip_disc = TClipDisc(sf_rna_algnmt, s_wfolder, n_jobs, sf_ref)
    # sf_disc_fa_tmp = self.s_wfolder + "temp_disc.fa"
    sf_clip_fa_tmp = s_wfolder + "temp_clip.fq"  ##

    ####collect clipped and disc reads
    sf_disc_fa = s_wfolder + "collected_disc.fa"
    bin_size = 50000000  # block size for parallelization
    xclip_disc.collect_clipped_disc_reads_of_given_list(l_record, bin_size, sf_clip_fa_tmp, sf_disc_fa)