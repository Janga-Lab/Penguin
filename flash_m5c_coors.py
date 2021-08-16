##################################################
import re
import os

def cigar_parse(cigar,start):
    insertions = []
    alignments = []
    aligns=[]

    pattern = re.compile('([MIDNSHPX=])')
    values = pattern.split(cigar)[:-1] ## turn cigar into tuple of values
    paired = (values[n:n+2] for n in range(0,len(values),2)) ## pair values by twos
    i = 0 ## alignment coordinate index
    g = start ## genomic coordinate index
    gstop=0
    for pair in paired:
        l = int(pair[0]) ## length of CIGAR event
        t = pair[1] ## type of CIGAR event
        if t == 'M': ## if match, return consecutive coordinates
            alignments.append((g, g+i+l-i,(i, i + l))) ## (genomic offset, (alignment.start, alignment.end))
            aligns.append((g, g+i+l-i, i, i+l)) ## (genomic offset, (alignment.start, alignment.end))
            i += l
            g += l
            gstop=g+i+l-i
        elif t == 'D': ## skip 'l' number of coordinates in reference
            g += l
        elif t == 'I': ## insertion of 'l' length
            insertions.append((i, i + l))
            i += l
        elif t == 'N': ## skipped region from the reference
            g += l
        elif t == 'S': ## soft clipping (clipped sequences present in SEQ)
            i += l
        elif t == 'H': ## hard clipping (clipped sequences NOT present in SEQ)
            pass
        elif t == 'P': ## padding (silent deletion from padded reference)
            pass
        elif t == '=': ## sequence match
            pass
        elif t == 'X': ## sequence mismatch
            pass
    return aligns, gstop


filename = input("Input the bed filename with its extension:") 
if filename: 
    locs=open(filename,'r')
else:
    print('Empty Bed file name')

sam_file=[]
'''
file3=open('hek.sam','r') #### Provide input sam file
for u in file3:
    u1=u.split("\t")
    if len(u1) > 3:
        sam_file.append(u1[0]+' '+u1[1]+' '+u1[2]+' '+u1[3]+' '+u1[4]+' '+u1[5]+' '+u1[6])

'''
file3=open('hek.sam','r') #### Provide input sam file
#file3=open('Hela_align_1.sam','r') #### Provide input sam file
for u in file3:
    #print(u)
    u1=u.split( )
    #print(len(u1))
    if len(u1)>3:
        sam_file.append(u1[0]+' '+u1[1]+' '+u1[2]+' '+u1[3]+' '+u1[4]+' '+u1[5]+' '+u1[6] )

total=len(sam_file)
print(total)
print("Sam file read successfully..starting coordinate extraction...")

mod_d=dict()

dict_len=0
#newopen = open('new_x.txt', 'w')
mods=[]
for e in locs:
    e1=e.split( )
    #print(e1[5])
    #newopen.write(e1[5]+'\n')
    if e1[5] == '+':  #change 2 to 5 as the sign loctes in 6th column 
        #print(e1[0])
        mods.append(e1[0]+'_'+e1[1])
        chr=''.join(e1[0])
        chr_loc=''.join(e1[1])
        if chr in mod_d:
            mod_d[chr].append(int(chr_loc))
        else:
            mod_d[chr] = [int(chr_loc)]
#newopen.close()
for key, value in mod_d.items() :
    # print (key, len(value))
    # print('\n')
    dict_len=dict_len+len(value)
    
#print(dict_len)
file2=open('Ps_Modification_coors_hek_complete.txt','a') ### define output filename
c=0
### loop through samfile ####
for t, i in enumerate(sam_file):
    if not i.startswith('@'):
        i1=i.split()
        s_chr=''.join(i1[2])
        #print(s_chr)
        s_chr_s=''.join(i1[3])
    #  print(i1[0]+' '+i1[3]+' '+i1[5])   
        coordinates,gstop=cigar_parse(i1[5],int(i1[3]))
        #print(i1[5])
        # if 'chr'+i1[2] == e1[0]:
        #print(coordinates)    
        # print ("chr"+i1[2]+' '+i1[3])
        #print(mod_d[s_chr])  #problem here what this should return?
        try:
            for v in mod_d[s_chr]:
                # print(v)
                # if v in range(int(i1[3]),int(i1[3])+gstop):
                if int(v)-2 >= int(i1[3]) and int(v)+2 <= int(i1[3])+gstop:
                    # print(s_chr_s+' '+str(v)+' '+str(int(i1[3])+gstop))
                    # print(gstop)
                    for r in coordinates:
                        if int(r[0]) <= int(v) and int(r[1]) >= int(v): 
                            # print(i1[3])
                            # print(e1[1])
                            # print(r) 
                            genomic_coor=list(range(int(r[0]),int(r[1])+1))
                            # print(genomic_coor)
                            seq_coor=list(range(int(r[2]),int(r[3])+1))
                            # print(seq_coor)
                            i_gc=genomic_coor.index(int(v))
                            i_sc=seq_coor[i_gc]
                            # print(str(s_chr)+' '+str(v)+' '+str(i_sc)+' '+str(i1[0]))
                            file2.write(str(s_chr)+' '+str(v)+' '+str(i_sc)+' '+str(i1[0])+'\n')
        except KeyError:
            continue

                    # file2.write(i+'\n')   
                    # print(e1[1]+' '+str(i1[3]))
    c=c+1
   
