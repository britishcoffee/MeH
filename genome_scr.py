
import random
import math
import pysam
import csv
import sys
import pandas as pd
import numpy as np
import datetime
import time as t


#---------------------------------------

# Functions definition

#---------------------------------------

# Check whether a window has enough reads for complete/impute
def enough_reads(window,w,complete):
    temp=np.isnan(window).sum(axis=1)==0
    if complete: # For heterogeneity estimation
        return temp.sum()>=2**(w-2)
    else:  # for imputation
        tempw1=np.isnan(window).sum(axis=1)==1
        return temp.sum()>=2**(w-1) and tempw1.sum()>0
    
def impute(window,w):
    full_ind=np.where(np.isnan(window).sum(axis=1)==0)[0]
    part_ind=np.where(np.isnan(window).sum(axis=1)==1)[0]
    for i in range(len(part_ind)):
        sam = []
        # which column is nan
        pos=np.where(np.isnan(window[part_ind[i],:]))[0]
        if np.unique(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos]).shape[0]==1:
            window[part_ind[i],pos]=window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos][0]
        else:
            #print("win_part i pos =",window[part_ind[i],pos])
            for j in range(len(full_ind)):
                if (window[part_ind[i],:]==window[full_ind[j],:]).sum()==w-1:
                    sam.append(j)
            if len(sam)>0:
                s1=random.sample(sam, 1)
                s=window[full_ind[s1],pos]
            else:
                s=random.sample(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos].tolist(), k=1)[0]
            window[part_ind[i],pos]=np.float64(s)
            #print("win_part i =",window[part_ind[i],pos])
            #print("s = ",np.float64(s))
    return window 
   
def getcomplete(window,w):
    temp=np.isnan(window).sum(axis=1)==0
    mat=window[np.where(temp)[0],:]
    return mat

def PattoDis(mat):
    s=mat.shape[0]
    dis=np.zeros((s,s))
    for i in range(s):
        for j in range(s):
            if j<i:
                dis[i,j]=dis[j,i]=Ham_d(mat.iloc[i,],mat.iloc[j,]) 
    return dis
        
def Ham_d(pat1,pat2): 
    return (pat1!=pat2).sum()

# input a window of w CGs and output a list of proportions with starting genomic location and genomic distance across
def window_summ(pat,start,dis,chrom): 
    m=np.shape(pat)[0]
    d=np.shape(pat)[1]
    all_pos=np.zeros((2**d,d))
    for i in range(d): 
        all_pos[:,i]=np.linspace(0,2**d-1,2**d)%(2**(i+1))//(2**i)
    #print(all_pos)
    prob=np.zeros((2**d,1))
    #print(prob)
    for i in range(2**d): 
        count = 0
        for j in range(m):
            if (all_pos[i,:]==pat.iloc[j,:]).sum()==d:
                count+=1
                #print(count)
        prob[i]=count
    if d==3:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p1':prob[0],'p2':prob[1],'p3':prob[2],'p4':prob[3],\
                    'p5':prob[4],'p6':prob[5],'p7':prob[6],'p8':prob[7],'dis':dis})    
    if d==4:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p1':prob[0],'p2':prob[1],'p3':prob[2],'p4':prob[3],\
                    'p5':prob[4],'p6':prob[5],'p7':prob[6],'p8':prob[7],'p9':prob[8],'p10':prob[9],\
                    'p11':prob[10],'p12':prob[11],'p13':prob[12],'p14':prob[13],'p15':prob[14],\
                    'p16':prob[15],'dis':dis})   
    if d==5:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p1':prob[0],'p2':prob[1],'p3':prob[2],'p4':prob[3],\
                    'p5':prob[4],'p6':prob[5],'p7':prob[6],'p8':prob[7],'p9':prob[8],'p10':prob[9],\
                    'p11':prob[10],'p12':prob[11],'p13':prob[12],'p14':prob[13],'p15':prob[14],\
                    'p16':prob[15],'p17':prob[16],'p18':prob[17],'p19':prob[18],'p20':prob[19],\
                    'p21':prob[20],'p22':prob[21],'p23':prob[22],'p24':prob[23],'p25':prob[24],\
                    'p26':prob[25],'p27':prob[26],'p28':prob[27],'p29':prob[28],'p30':prob[29],\
                    'p31':prob[30],'p32':prob[31],'dis':dis})
    return out

def MeHperwindow(pat,start,dis,chrom,D,all_pos,w,q=2): 
    prob=np.zeros((2**w,1))
    m=np.shape(pat)[0]
    #print(prob)
    for i in range(2**w): 
        count = 0
        for j in range(m):
            if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                count+=1
        prob[i]=count
    interaction=np.multiply.outer(prob/m,prob/m).reshape((2**w,2**w))
    Q=sum(sum(D*interaction))
    if Q==0:
        div=0
    else:
        div=(sum(sum(D*(interaction**2)))/Q)**(-0.5)
    
    out=pd.DataFrame({'chrom':chrom,'pos':start,'MeH':div,'dis':dis}, index=[0])    
    return out
    
def PW_H(pat,q=2,w=3,dist=Ham_d,prop=False): 
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    prob=np.zeros((2**w,1))
    #print(prob)
    if not prop:
        m=np.shape(pat)[0]
        for i in range(2**w): 
            count = 0
            for j in range(m):
                if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                    count+=1
                    #print(count)
            prob[i]=count
    if prop:
        prob=pat
    
    D=PattoDis(pd.DataFrame(all_pos))
    m=prob.sum(axis=0)
    interaction=np.multiply.outer(prob/m,prob/m).reshape((2**w,2**w))
    Q=sum(sum(D*interaction))
    if Q==0:
        return Q
    else:
        div=sum(sum((D*(interaction/Q)**2)))**(1/(2*(1-q)))
        return div
    
def Ab_H(pat,q=2,w=3,dist=Ham_d,prop=False): 
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    prob=np.zeros((2**w,1))
    #print(prob)
    if not prop:
        m=np.shape(pat)[0]
        for i in range(2**w): 
            count = 0
            for j in range(m):
                if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                    count+=1
                    #print(count)
            prob[i]=count     
    if prop:
        prob=pat
    m=prob.sum(axis=0)
    out=(((prob/m)**2).sum(axis=0))**(1/(2*(1-q)))
    return np.float64(out)

def Ent_H(pat,q=2,w=3,dist=Ham_d,prop=False): 
    if not prop:
        all_pos=np.zeros((2**w,w))
        for i in range(w): 
            all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
        prob=np.zeros((2**w,1))
        #print(prob)
        m=np.shape(pat)[0]
        for i in range(2**w): 
            count = 0
            for j in range(m):
                if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                    count+=1
                    #print(count)
            prob[i]=count/m
    if prop:
        prob=pat
        m=prob.sum(axis=0)
        prob=prob/m
    prob=prob.to_numpy()
    out=0
    for i in prob:
        if i>0:
            out-=i*math.log2(i)/w
    return np.float64(out)

def Epi_H(pat,q=2,w=3,dist=Ham_d,prop=False): 
    if not prop:
        all_pos=np.zeros((2**w,w))
        for i in range(w): 
            all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
        prob=np.zeros((2**w,1))
        #print(prob)
        #print(all_pos)
        m=np.shape(pat)[0]
        
        for i in range(2**w): 
            count = 0
            for j in range(m):
                if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                    count+=1
                    #print(count)
            prob[i]=count
    if prop:
        prob=pat
    m=prob.sum(axis=0)
    out=1-((prob/m)**2).sum(axis=0)
    #print(type(np.float64(out)))
    return np.float64(out)


def Phy_H(pat,q=2,dist="hamming_d",w=4,prop=False):
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    prob=np.zeros((2**w,1))
    #print(prob)
    if not prop:
        m=np.shape(pat)[0]
        for i in range(2**w): 
            count = 0
            for j in range(m):
                if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
                    count+=1
                    #print(count)
            prob[i]=count
    if prop:
        prob=np.append([0],pat)
    if dist=="hamming_d" and w==4:
        phylotree=np.append(np.append(np.append(np.append([0],np.repeat(0.5,16)),np.repeat(0.25,6)),[0.5]),np.repeat(0.25,6))
        #phylotree=np.repeat(0,1).append(np.repeat(0.5,16)).append(np.repeat(0.25,6)).append(0.5).append(np.repeat(0.25,6))
        count=np.zeros(30)
        #count<-rep(0,29)
        count[1:17]=prob[1,9,5,3,2,13,11,10,7,6,4,15,14,12,8,16]*2**w
        count[17]=count[4]+count[7]
        count[18]=count[9]+count[12]
        count[19]=count[1]+count[2]
        count[20]=count[3]+count[6]
        count[21]=count[17]+count[18]
        count[22]=count[19]+count[20]
        count[23]=count[21]+count[22]
        count[24]=count[5]+count[8]
        count[25]=count[10]+count[13]
        count[26]=count[24]+count[25]
        count[27]=count[23]+count[26]
        count[28]=count[11]+count[14]
        count[29]=count[27]+count[28]
        #Q=sum(sum(phylotree*count))    
    if dist=="WDK_d" and w==4: 
        phylotree=np.append(np.append(np.append(np.append([0],np.repeat(3,16)),np.repeat(1.5,6)),[3.2,0.8]),np.repeat(2,3),np.repeat(1.5,2))
        #phylotree=c(rep(3,16),rep(1.5,6),3.2,0.8,rep(2,3),1.5,1.5)
        count=np.zeros(30)
        #print(prob)
        count[1:17]=prob[1,9,5,3,2,13,11,10,7,6,4,15,14,12,8,16]
        count[17]=count[1]+count[2]
        count[18]=count[5]+count[8]
        count[19]=count[3]+count[6]
        count[20]=count[10]+count[13]
        count[21]=count[4]+count[7]
        count[22]=count[11]+count[14]
        count[23]=count[17]+count[18]
        count[24]=count[21]+count[22]
        count[25]=count[19]+count[20]
        count[26]=count[23]+count[24]
        count[27]=count[25]+count[26]
        count[28]=count[9]+count[12]
        count[29]=count[27]+count[28]
        #Q=sum(phylotree*count)
    if dist=="WDK_d" and w==3:
        phylotree=np.append(np.append(np.append([0],np.repeat(1.5,8)),np.repeat(0.75,3)),np.repeat(1.5,0.75))
        #phylotree=np.array(0).append(np.repeat(1.5,8)).append(np.repeat(0.75,3)).append(1.5,0.75)
        #phylotree=c(rep(1.5,8),rep(0.75,3),1.5,0.75)
        count=np.zeros(14)
        count[1:9]=prob[1:9]
        count[9]=count[1]+count[2]
        count[10]=count[5]+count[6]
        count[11]=count[3]+count[4]
        count[12]=count[9]+count[10]
        count[13]=count[11]+count[12]
        #Q=sum(phylotree*count)
    if dist=="hamming_d" and w==3:
        phylotree=np.append(np.append(np.append([0],np.repeat(0.5,8)),np.repeat(0.25,3)),[0.5,0.25])
        #phylotree=np.array(0).append(np.repeat(0.5,8)).append(np.repeat(0.25,3)).append(0.5,0.25)
        count=np.zeros(14)
        count[1:9]=prob[1:9]
        count[9]=count[1]+count[2]
        count[10]=count[5]+count[6]
        count[11]=count[3]+count[4]
        count[12]=count[9]+count[10]
        count[13]=count[11]+count[12]
        #print("count = ",count)
        #print("phylotree = ",phylotree)
    Q=sum(phylotree*count)
    div=sum(phylotree*((count/Q)**q))**(1/(1-q))
    return div

# Function : Result -> MeH
def MeH(Result,prop=True,fun=PW_H,w=3):
    out=pd.DataFrame(columns=['pos','MeH','dis'])
    for i in range(Result.shape[0]):
        a=fun(pat=pd.DataFrame(Result.iloc[i,1:1+2**w]),prop=prop)
        #print(a)
        toappend=pd.DataFrame(data={'MeH':a,'pos': Result.iloc[i,0], 'dis': Result.iloc[i,1+2**w]}, index=[i])
        out=out.append(toappend)
    return out

def impute(window,w):
    full_ind=np.where(np.isnan(window).sum(axis=1)==0)[0]
    part_ind=np.where(np.isnan(window).sum(axis=1)==1)[0]
    for i in range(len(part_ind)):
        sam = []
        # which column is nan
        pos=np.where(np.isnan(window[part_ind[i],:]))[0]
        if np.unique(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos]).shape[0]==1:
            window[part_ind[i],pos]=window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos][0]
        else:
            #print("win_part i pos =",window[part_ind[i],pos])
            for j in range(len(full_ind)):
                if (window[part_ind[i],:]==window[full_ind[j],:]).sum()==w-1:
                    sam.append(j)
            if len(sam)>0:
                s1=random.sample(sam, 1)
                s=window[full_ind[s1],pos]
            else:
                s=random.sample(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos].tolist(), k=1)[0]
            window[part_ind[i],pos]=np.float64(s)
            #print("win_part i =",window[part_ind[i],pos])
            #print("s = ",np.float64(s))
    return window 


def genome_scr(bamfile,w,fa,silence=False):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    aggreC = pd.DataFrame(columns=['Qname'])
    Result1 = pd.DataFrame(columns=['chrom','pos','MeH','dis'])
    never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    D=PattoDis(pd.DataFrame(all_pos))
    start=datetime.datetime.now()
    for pileupcolumn in samfile.pileup():
        chrom = pileupcolumn.reference_name
        if not silence: :
            if (pileupcolumn.pos % 50000 < 1):
                print("%s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,Result1.shape[0]))
        if (fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+2)=='CG'):        
            temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos])
            pileupcolumn.set_min_base_quality(0)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                    d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                    df2 = pd.DataFrame(data=d)
                    #df2.head()
                    temp=temp.append(df2, ignore_index=True)
            #temp.head()
            if (not temp.empty):
                #temp.head()
                aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                aggreC = aggreC.drop_duplicates()
            
        if Result1.shape[0]==0 and aggreC.shape[1] == (2*w):
            never = False
            aggreC = aggreC.replace(['C','G'],1)
            aggreC = aggreC.replace(['A','T'],0)
            aggreC = aggreC.replace(['N'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                if i<w and enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                    dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                    chrom=chrom,D=D,all_pos=all_pos,w=w)
                    Result1=Result1.append(toappend)
            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
        #------------------
        #  SECONDARY CASE
        #------------------
        if (aggreC.shape[1] == (3*w-1)):
            aggreC = aggreC.replace(['C','G'],1)
            aggreC = aggreC.replace(['A','T'],0)
            aggreC = aggreC.replace(['N'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                if i>w-2 and i<2*w and enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                chrom=chrom,D=D,all_pos=all_pos,w=w)
                    Result1=Result1.append(toappend)
                    if Result1.shape[0] % 5000 == 0:   
                        Result1.to_csv(r"MeHdata/gs_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,Result1.shape[0],pileupcolumn.pos))
            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)      
    if Result1.shape[0]>0:   
        Result1.to_csv(r"MeHdata/gs_%s.csv"%(filename),index = False, header=True)
        print("Done. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,Result1.shape[0],pileupcolumn.pos))
    
def split_bam(samplenames,Folder): 
        # get bam size
        spbam_list = []
        bamfile = samplenames + '.bam'
        statinfo_out = os.stat(Folder+bamfile)
        bamsize = statinfo_out.st_size
        samfile = pysam.Samfile(Folder+bamfile, "rb")
        fileout_base = os.path.splitext(bamfile)[0] # filename
        ext = '.bam'
        x = 0
        fileout = Folder+fileout_base+"_" + str(x)+ext # filename_x.bam
        print("fileout ",fileout)
        header = samfile.header
        
        outfile = pysam.Samfile(fileout, "wb", header = header)
        sum_Outfile_Size=0
        for reads in samfile.fetch():
            # ordered?
            outfile.write(reads)
            statinfo_out = os.stat(fileout)
            outfile_Size = statinfo_out.st_size
            #if(outfile_Size >=1073741824 and sum_Outfile_Size <= infile_Size):
            if(outfile_Size >=57374182 and sum_Outfile_Size <= bamsize):
                sum_Outfile_Size = sum_Outfile_Size + outfile_Size
                x = x + 1
                spbam_list.append(fileout_base + "_" + str(x)+ext)
                outfile.close()
                pysam.index(fileout)
                fileout = Folder+fileout_base + "_" + str(x)+ext
                print("fileout ",fileout)
                outfile = pysam.Samfile(fileout, "wb",header = header)
            
        outfile.close()
        pysam.index(fileout)
        
import argparse
parser = argparse.ArgumentParser()
#parser.add_argument("-i", "--input", dest="filename", type=validate_file, required=True, help="enter input file", metavar="FILE")
#parser.add_argument('-o', metavar='out-file', type=argparse.FileType('wt'))
#parser.add_argument("-f", "--filename",type=str,help='input file prefix',required=True)
parser.add_argument("-w", "--windowsize",type=int,help='number of CGs')
#parser.add_argument("-c", "--chromosome",type=str,help='chromosome')
parser.add_argument("-c", "--cores",type=str,help='number of cores')
args = parser.parse_args()

import sys
import os
import pandas as pd
import multiprocessing
from joblib import Parallel, delayed

#num_cores = multiprocessing.cpu_count()
                                                
if __name__ == "__main__":

    Folder = 'MeHdata/'

    files = os.listdir(Folder)
    bam_list = []
    # all samples' bam files
    for file in files: 
        filename, file_extension = os.path.splitext(file)
        if file_extension == '.fa':
            fa = filename
        if file_extension == '.bam':
            bam_list.append(filename)
            
    if 'cores' in args: 
        num_cores = args.cores
    else:
        num_cores = 4
        
    Parallel(n_jobs=num_cores)(delayed(split_bam)(bamfile,Folder=Folder) for bamfile in bam_list)
    
    spbam_list = []
    tempfiles = os.listdir(Folder)
    for file in tempfiles:
        filename, file_extension = os.path.splitext(file)
        if file_extension=='.bam' and filename not in bam_list:
            spbam_list.append(filename)
    print(spbam_list)
        
    start=t.time()
    Parallel(n_jobs=args.cores)(delayed(genome_scr)(bamfile,w=args.windowsize,fa=fa) for bamfile in spbam_list)
    print(t.time()-start)

    # merge .csv within sample
    for file in spbam_list:
        #filename, file_extension = os.path.splitext(file)
        sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + str(sample) + '.csv'
                toapp_dir = Folder + 'gs_' +file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    #os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend.to_csv(res_dir,index = False,header=True)
                    #os.remove(toapp_dir)

    #os.chdir('../')
    #os.chdir(outputFolder)
    
    # merge .csv between samples
    for sample in bam_list: 
        tomerge_dir = Folder + str(sample) + '.csv' 
        res_dir = Folder + 'Results.csv'
        if os.path.exists(res_dir):
            Result = pd.read_csv(res_dir)
            Tomerge = pd.read_csv(tomerge_dir)
            Tomerge = Tomerge.drop(columns='dis')
            Tomerge = Tomerge.rename(columns={'MeH': sample})
            Result = Result.merge(Tomerge, on=['chrom','pos'])
            Result = Result.drop_duplicates() 
            Result.to_csv(Folder+'Results.csv',index = False,header=True)
            os.remove(tomerge_dir)
        else:
            Result = pd.read_csv(tomerge_dir)
            Result = Result.drop(columns='dis')
            Result = Result.rename(columns={'MeH': sample})
            Result.to_csv(Folder+'Results.csv',index = False,header=True)
            os.remove(tomerge_dir)

    Result.to_csv(Folder+'Results.csv' ,index = False,header=True)

    for filename in spbam_list:
        file = Folder + filename + '.bam'
        os.remove(Folder + file)


    
