'''
Created on Sep 29, 2011

@author: ksahlin
'''

class scaffold:
    '''
    classdocs
    '''


    def __init__(self, scaffold_name, scaffold_contigs,scaffold_length, scaffold_left_nbrs, scaffold_right_nbrs ):
        '''
        Constructor
        '''
        self.name = scaffold_name           #String name
        self.contigs = scaffold_contigs     #list of contig objects that are ordered ass they should be placed in the scaffold
        self.s_length = scaffold_length     #integer of total length of the scaffold   
        self.left_nbrs_obs = scaffold_left_nbrs
        self.right_nbrs_obs = scaffold_right_nbrs
        
        #self.stop= len(scaffold_contigs)
        #self.index=0
        #self.positions = scaffold_positions
        #self.directions = scaffold_directions
        
                
#    #standard forward iterator over contig objects to change positions directions
#    def __iter__(self):
#        return self
#    def next(self):
#        if self.index == self.stop:
#            raise StopIteration
#        self.index = self.index + 1
#        return self.contigs[self.index-1]
#        
#    #Reverse iterator over contig objects to change positions directions
#    def reviter(self):
#        return self
#    def prev(self):
#        if self.index == 0:
#            raise StopIteration
#        self.index = self.index - 1
#        return self.contigs[self.index]
#    