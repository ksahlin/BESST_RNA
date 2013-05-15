'''
Created on Oct 4, 2011

@author: ksahlin
'''
def PrintOutHaplotypes(Haplotypes,Contigs,output_dest,C_dict):
    haplotype_file=open(output_dest+'/haplotypes.fa','a')
    for cont_hapl in Haplotypes: 
        print >>haplotype_file, '>'+ Haplotypes[cont_hapl][0]+ ' (variant of '+cont_hapl+')'
        hapl_contig = C_dict[Haplotypes[cont_hapl][0]]
        for i in range(0,len(hapl_contig),60):
            print >>haplotype_file, hapl_contig[i:i+60]
        #remove the outputted haplotype, the other one is left in graph
        del Contigs[Haplotypes[cont_hapl][0]]
    return(Contigs)

def PrintOutRepeats(Repeats,Contigs,output_dest,C_dict):
    repeat_file=open(output_dest+'/repeats.fa','w')
    for cont_obj in Repeats: 
        print >>repeat_file, '>'+cont_obj.name 
        repeat_contig = C_dict[cont_obj.name]
        for i in range(0,len(repeat_contig),60):
            print >>repeat_file, repeat_contig[i:i+60]
        del Contigs[cont_obj.name]
    return(Contigs)

    
def WriteToF(F,Contigs,list_of_contigs):
    info_list=[]
    for cont_obj in list_of_contigs:
        info_list.append((cont_obj.name, cont_obj.direction,cont_obj.position, cont_obj.length,cont_obj.links))
        if cont_obj.position < 0:
            print 'Write to F: Position is negative!', cont_obj.position, cont_obj.name, cont_obj.direction
        del Contigs[cont_obj.name]        
    F.append(info_list)
    return(Contigs,F)

def PrintOutput(F,C_dict,Information,output_dest,param,pass_nr):
    import os
    try:
        os.mkdir(param.output_directory +'/pass'+str(pass_nr))
    except OSError:
        #directory is already created
        pass
    contigs_before=len(C_dict)
    contigs_after=len(F)
    out_scaf_file=open(param.output_directory +'/pass'+str(pass_nr)+'/New_RNA_Scaffolds-pass'+str(pass_nr),'w')
    print >>Information, 'Contigs before scaffolding: '+ str(contigs_before)+'\n'
    print >>Information, '(super)Contigs after scaffolding: '+ str(contigs_after)+'\n'
    gff_file=open(param.output_directory +'/pass'+str(pass_nr)+'/info-pass'+str(pass_nr)+'.gff','w')
    AGP_file=open(param.output_directory +'/pass'+str(pass_nr)+'/info-pass'+str(pass_nr)+'.agp','w')
    print >>gff_file, '#gff-version 3'
    print >>AGP_file, '#APG file\n#lw-scaffolder output'    
    output=open(param.output_directory +'/pass'+str(pass_nr)+'/Scaffolds-pass'+str(pass_nr)+'.fa','w')
    header_index=0
    prev_pos=0
    print_list=[]
    print >>out_scaf_file, 'New scaffold name\tincluded old scaffolds\n '
    component_count=0
    for scaf_ in reversed(F):
        #sort contigs in scaf w.r.t position here
        scaf=sorted(scaf_, key=lambda tuple: tuple[2])
        header_index+=1
        print_list=[]
        scaf_len=len(scaf)
        if scaf_len > 1:
            scaf_name = 'new_scaffold_'+str(header_index)
            print >>output, '>'+scaf_name
        #print_list.append('>scaffold_'+str(header_index)+'\n')
        else:
            scaf_name = scaf[0][0]
            print >>output, '>'+scaf_name
        
        prev_pos=0
        component_count=0
        for i in range (0,scaf_len):
            name=scaf[i][0]
            direction=scaf[i][1]
            pos=scaf[i][2]
            length=scaf[i][3]
            
            #Get the number of left and right links            
            if i-1>=0:
                left_name=scaf[i-1][0]
                try:
                    nr_links_left=scaf[i][4][left_name]
                except KeyError:
                    nr_links_left='!'                    
            else:
                nr_links_left=0
                
            if i+1<=scaf_len-1:
                right_name=scaf[i+1][0]
                try:
                    nr_links_right=scaf[i][4][right_name]
                except KeyError:
                    nr_links_right='!'                    
            else:
                nr_links_right=0
                
            gap=pos-prev_pos
            N_seq='n'*gap
            if gap < 0: #remove a part of the last contig so we can append the next one (overwrite it's sequence)
                last_contig=print_list.pop()
                last_contig=last_contig[0:gap]
                print_list.append(last_contig)
                
            print_list.append(N_seq)            
            if direction: 
                print_list.append(C_dict[name])                          
            else:  #needs to be reverse complemented before outputted
                rev_comp=RevComp(C_dict[name])
                print_list.append(rev_comp)  
                         
            sign='+' if direction else '-'
            #FOR GFF file
            if scaf_len > 1 and  i>0 and gap > 0:
                print >>gff_file, scaf_name +'\tBESST_RNA \tfragment\t' + str(prev_pos+1)+'\t'+str(pos)+'\t.\t'+sign+'\t.\t'+'Parent='+scaf_name
            #FOR AGP file            
            if i>0 and gap > 0: 
                component_count+=1
                print >>AGP_file, scaf_name + '\t'+str(prev_pos+1)+'\t'+str(pos)+'\t'+str(component_count)+'\t'+'N\t'+str(gap)+'\tfragment\tyes\t'            
            component_count+=1            
            prev_pos=pos+length  
            print >>AGP_file, scaf_name + '\t'+str(pos+1)+'\t'+str(prev_pos)+'\t'+str(component_count)+'\t'+'W\t'+name+'\t1\t'+str(length)+'\t'+sign  
            
            #FOR GFF file
            if scaf_len > 1:              
                print >>gff_file, scaf_name + '\tBESST_RNA\tcontig\t' +str(pos+1)+'\t'+str(prev_pos)+'\t.\t'+sign +'\t.\t'+'Parent='+scaf_name+';ID='+name +';Left-links='+str(nr_links_left) +';Right-links='+str(nr_links_right) 

        if scaf_len > 1:
            print >>out_scaf_file, scaf_name + '\t' +''.join([obj[0]+', ' for obj in scaf])[:-2] 
            print >>gff_file, scaf_name + '\tBESST_RNA\tscaffold\t' +'1'+'\t'+str(prev_pos)+'\t.\t'+ '+' +'\t.\t'+'ID='+scaf_name              
        items=len(print_list)                
        output_scaf=''.join([print_list[i] for i in range(0,items)])
        for i in range(0,len(output_scaf),60):
            print >>output, output_scaf[i:i+60]

        #print >>output, ''.join([print_list[i] for i in range(0,items)])
    return()


def RevComp(string):
    #rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X'}
    ## this is the madness reverse complement table that needs to be specified when working with AbySS...
    rev_nuc={'A':'T','C':'G','G':'C','T':'A','a':'t','c':'g','g':'c','t':'a','N':'N','X':'X','n':'n','Y':'R','R':'Y','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','H':'D','D':'H','y':'r','r':'y','k':'m','m':'k','s':'s','w':'w','b':'v','v':'b','h':'d','d':'h'}
    rev_comp=''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)
            

