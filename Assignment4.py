#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 15:44:54 2022

@author: kurup
"""
import numpy as np
from math import factorial

###################################### if you have an explicit folder for files, mention its navigation in the variable "path"
path="" # <----- give the explicit location to all supporting files
##############################################################################################################################
#                   INITIALISATION OF VARIABLES AND READING THE FILES

file = open(path+"chrX_map.txt",'r')
map = file.readlines()
map=np.array([int(a.split("\n")[0]) for a in map])
file.close()

file = open(path+"chrX.fa",'r')
file.readline()
dot_fa = file.readlines()
file.close()

file = open(path+"chrX_last_col.txt",'r')
last_col = file.readlines()
file.close()
A=np.zeros(len(last_col))
C=A.copy()
G=A.copy()
T=A.copy()
A[0],C[0],G[0],T[0]=last_col[0].count("A"),last_col[0].count("C"),last_col[0].count("G"),last_col[0].count("T")
for i in range(1,len(last_col)):
    A[i],C[i],G[i],T[i]=last_col[i].count("A")+A[i-1],last_col[i].count("C")+C[i-1],last_col[i].count("G")+G[i-1],last_col[i].count("T")+T[i-1]

file = open(path+'reads','r')
file.readline()
reads = file.readlines()
file.close()
#  sbatch submit_job.sh
RXon = np.zeros(6)
GXon = np.zeros(6)
Xon = np.zeros(6)
RX_ind=np.array([[149249757, 149249868],[149256127, 149256423],[149258412, 149258580],[149260048, 149260213],[149261768, 149262007],[149264290, 149264400]])
GX_ind=np.array([[149288166, 149288277],[149293258, 149293554],[149295542, 149295710],[149297178, 149297343],[149298898, 149299137],[149301420, 149301530]])
#######################################################################################################
def search_support(lis,start,start_index,end,end_index,character): # PROCEDURE SUPPORTING SS
    s_1 = lis[start]-last_col[start][start_index:].count(character)+1 #
    s_2 = lis[end]-last_col[end].count(character)+last_col[end][:end_index+1].count(character)
    return s_1,s_2
    
def ss(inp):  #SEARCH SEQUENCE PROCEDURE
    ret = [0,A[-1]+C[-1]+G[-1]+T[-1]+1]
    seq=list(range(len(inp)))
    seq.reverse()

    for i in seq:        
        character = inp[i]
        start = int(ret[0]//100)
        end = int(ret[1]//100)
        start_index=int(ret[0]%100)
        end_index=int(ret[1]%100)
        if character=='A':
            s_1,s_2=search_support(A,start,start_index,end,end_index,character)
            if s_1>s_2:
                return ret,i+1
            ret = [s_1-1,s_2-1]
        elif character=='C':
            s_1,s_2=search_support(C,start,start_index,end,end_index,character)
            if s_1>s_2:
                return ret,i+1
            ret = [A[-1]+s_1-1,A[-1]+s_2-1]
        elif character=='G':
            s_1,s_2=search_support(G,start,start_index,end,end_index,character)
            if s_1>s_2:
                return ret,i+1
            ret = [A[-1]+C[-1]+s_1-1,A[-1]+C[-1]+s_2-1]
        elif character=='T':
            s_1,s_2=search_support(T,start,start_index,end,end_index,character)
            if s_1>s_2:
                return ret,i+1            
            ret = [A[-1]+C[-1]+G[-1]+s_1-1, A[-1]+C[-1]+G[-1]+s_2-1]
    return ret,0

#######################################################################################################
def sub_step(read,R,G):
    ss_ret,move_to = ss(read)
    for i in range(int(ss_ret[0]),int(ss_ret[1])+1):
        id = int(map[i])-move_to           
        ref=""
        len_read=len(read)
        block = int(id//100)
        offset = int(id%100)
        ref += dot_fa[block][offset:-1]
        while len(ref)<len_read:
            block+=1
            ref += dot_fa[block][:-1]
        ref = ref[:len_read]

        not_matched = 0
        allowed_buffer_match=2
        if len(read)==len(ref):
            for i in range(len(read)):
                if read[i]!=ref[i]:
                    not_matched+=1
                    if not_matched>allowed_buffer_match:
                        return R,G
            
            for ex_ind in range(6):
                if RX_ind[ex_ind][0]<= id <= RX_ind[ex_ind][1]:
                    R[ex_ind]=1
                if GX_ind[ex_ind][0]<= id <= GX_ind[ex_ind][1]:
                    G[ex_ind]=1
    return R,G

def step(read): # EXECUTED FOR EVERY LINE READ FROM THE reads FILE
    R=np.zeros(6)
    G=np.zeros(6)
    read = read[:-1]
    read=read.replace('N','A')
    comp_read=read[::-1]
    comp_read=comp_read.replace('A','1')
    comp_read=comp_read.replace('T','2')
    comp_read=comp_read.replace('C','3')
    comp_read=comp_read.replace('G','4')
    comp_read=comp_read.replace('1','T')
    comp_read=comp_read.replace('2','A')
    comp_read=comp_read.replace('3','G')
    comp_read=comp_read.replace('4','C')
   
    R,G=sub_step(read,R,G)
    R,G=sub_step(comp_read,R,G)
    for i in range(6):            
        if R[i]==G[i] and R[i]==1:
            Xon[i]+=0.5
            RXon[i]+=1
            GXon[i]+=1
        elif R[i]==1 or G[i]==1:
            Xon[i]+=1
            if R[i]==1:
                RXon[i]+=1
            if G[i]==1:
                GXon[i]+=1
for r in reads:
    step(r)
def comb(a,b):
    return factorial(int(a))/((factorial(int(b)))*(factorial(int(a-b))))

#                   RESULTS
#print(RXon)
#print(GXon)
ret_prob = np.zeros(4)
probabs=[np.array([(1/3),(1/3),(1/3),(1/3)]),
         np.array([(1/2),(1/2),(0.0),(0.0)]),
         np.array([(1/4),(1/4),(1/2),(1/2)]),
         np.array([(1/4),(1/4),(1/4),(1/2)])]
for iter_Exon in range(1,5):
    r=RXon[iter_Exon]
    g=GXon[iter_Exon]
    for iter_pro in range(4):
        p=probabs[iter_pro][iter_Exon-1]
        q=1-p
        add_on=comb((r+g),r)*(p**r)*(q**g)
        ret_prob[iter_pro] += float(add_on)
probs=ret_prob/ret_prob.sum()
print(f'Exon           : {Xon}')
print(f'Probabilities  : {probs}')
print(f'Max Probability: {probs.max()}')
print(f'Slot Number    : {int(np.where(probs==probs.max())[0])+1}')