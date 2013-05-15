'''
Created on Mar 31, 2012

@author: ksahlin
'''

class parameter(object):
    '''
    classdocs
    '''


    def __init__(self,parameter_mean_coverage = None, parameter_std_dev_coverage = None, parameter_mean_ins_size = None , parameter_std_dev_ins_size = None, 
                 parameter_output_directory = None, parameter_bamfile = None, parameter_read_len = None, parameter_rel_weight = None, parameter_ins_size_threshold = None
                 , parameter_contigfile = None, parameter_edgesupport = None, parameter_contig_threshold = None, parameter_scaffold_indexer = 0, parameter_first_lib = None,
                 parameter_cov_cutoff = None, parameter_tot_assembly_length = None, parameter_current_NG50 = None, parameter_current_LG50 = None,
                 parameter_hapl_ratio = None, parameter_hapl_threshold = None, parameter_detect_haplotype = None, parameter_detect_duplicate = None, parameter_gff_file = None,
                 parameter_information_file = None, parameter_fosmidpool = None):
        
        self.mean_coverage = parameter_mean_coverage
        self.std_dev_coverage = parameter_std_dev_coverage
        self.mean_ins_size = parameter_mean_ins_size
        self.std_dev_ins_size = parameter_std_dev_ins_size
        self.output_directory = parameter_output_directory
        self.bamfile = parameter_bamfile
        self.read_len = parameter_read_len
        self.rel_weight = parameter_rel_weight
        self.ins_size_threshold = parameter_ins_size_threshold
        self.contigfile = parameter_contigfile
        self.edgesupport = parameter_edgesupport
        self.contig_threshold = parameter_contig_threshold
        self.scaffold_indexer = parameter_scaffold_indexer
        self.first_lib = parameter_first_lib
        self.cov_cutoff = parameter_cov_cutoff
        self.tot_assembly_length = parameter_tot_assembly_length
        self.current_NG50 = parameter_current_NG50
        self.current_LG50 = parameter_current_LG50
        self.hapl_ratio = parameter_hapl_ratio
        self.hapl_threshold = parameter_hapl_threshold 
        self.detect_haplotype = parameter_detect_haplotype
        self.detect_duplicate = parameter_detect_duplicate
        self.gff_file = parameter_gff_file
        self.information_file = parameter_information_file
        self.fosmidpool = parameter_fosmidpool
        
        '''
        Constructor
        '''
        