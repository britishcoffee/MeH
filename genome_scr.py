##
#---------------------------------------------------------------------
# SERVER only input all files (.bam and .fa) output MeH matrix in .csv
# May 17, 2021
# FINAL github
#---------------------------------------------------------------------

import random
import math
import pysam
import csv
import sys
import pandas as pd
import numpy as np
import datetime
import time as t
from collections import Counter, defaultdict, OrderedDict


#---------------------------------------

# Functions definition

#---------------------------------------

# Count # of windows with enough reads for complete/impute
def coverage(methbin,complete,w):
    count=0
    tot = 0
    meth=methbin.iloc[:,methbin.columns!='Qname']
    if len(meth.columns)>=w:
        for i in range(len(meth.columns)-w+1):
            # extract a window
            temp = meth.iloc[:,i:i+w].copy()
            #print(temp)
            tot = tot+1
            if (enough_reads(window=temp,complete=complete,w=w)):
                count=count+1
                #toprint=temp.notnull().sum(axis=1)>=w
                #print(toprint.sum())
        #print(count)
        #print(tot)
        return count/tot*100
    else: 
        return 0

# Check whether a window has enough reads for complete/impute
def enough_reads(window,w,complete):
    temp=np.isnan(window).sum(axis=1)==0
    if complete: # For heterogeneity estimation
        return temp.sum()>=2**(w-2)
    else:  # for imputation
        tempw1=np.isnan(window).sum(axis=1)==1
        return temp.sum()>=2**(w-2) and tempw1.sum()>0
    

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
    #temp=window.notnull().sum(axis=1)>=w
    #mat=window.iloc[np.where(temp)[0],:]
    #else:
    #    temp=mat.notnull().sum(axis=1)>=w-1
    return mat

def PattoDis(mat,dist=1):
    s=mat.shape[0]
    dis=np.zeros((s,s))
    for i in range(s):
        for j in range(s):
            if j<i:
                if dist==1: 
                    d=Ham_d(mat.iloc[i,],mat.iloc[j,]) 
                else:
                    d=WDK_d(mat.iloc[i,],mat.iloc[j,]) 
                dis[i,j]=dis[j,i]=d
    return dis
        
def Ham_d(pat1,pat2): 
    return (pat1!=pat2).sum()

def WDK_d(pat1,pat2): 
    d=0
    w=pat1.shape[0]
    for i in range(w): # k-1
        for j in range(w-i): # starting pos
            s=(w-i-1)*(1-np.all(pat1[j:j+i+1]==pat2[j:j+i+1]))
            d+=s
    return d

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
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p01':prob[0],'p02':prob[1],'p03':prob[2],'p04':prob[3],\
                    'p05':prob[4],'p06':prob[5],'p07':prob[6],'p08':prob[7],'dis':dis})    
    if d==4:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p01':prob[0],'p02':prob[1],'p03':prob[2],'p04':prob[3],\
                    'p05':prob[4],'p06':prob[5],'p07':prob[6],'p08':prob[7],'p09':prob[8],'p10':prob[9],\
                    'p11':prob[10],'p12':prob[11],'p13':prob[12],'p14':prob[13],'p15':prob[14],\
                    'p16':prob[15],'dis':dis})   
    if d==5:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p01':prob[0],'p02':prob[1],'p03':prob[2],'p04':prob[3],\
                    'p05':prob[4],'p06':prob[5],'p07':prob[6],'p08':prob[7],'p09':prob[8],'p10':prob[9],\
                    'p11':prob[10],'p12':prob[11],'p13':prob[12],'p14':prob[13],'p15':prob[14],\
                    'p16':prob[15],'p17':prob[16],'p18':prob[17],'p19':prob[18],'p20':prob[19],\
                    'p21':prob[20],'p22':prob[21],'p23':prob[22],'p24':prob[23],'p25':prob[24],\
                    'p26':prob[25],'p27':prob[26],'p28':prob[27],'p29':prob[28],'p30':prob[29],\
                    'p31':prob[30],'p32':prob[31],'dis':dis})
    if d==6:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'p01':prob[0],'p02':prob[1],'p03':prob[2],'p04':prob[3],\
                    'p05':prob[4],'p06':prob[5],'p07':prob[6],'p08':prob[7],'p09':prob[8],'p10':prob[9],\
                    'p11':prob[10],'p12':prob[11],'p13':prob[12],'p14':prob[13],'p15':prob[14],\
                    'p16':prob[15],'p17':prob[16],'p18':prob[17],'p19':prob[18],'p20':prob[19],\
                    'p21':prob[20],'p22':prob[21],'p23':prob[22],'p24':prob[23],'p25':prob[24],\
                    'p26':prob[25],'p27':prob[26],'p28':prob[27],'p29':prob[28],'p30':prob[29],\
                    'p31':prob[30],'p32':prob[31],'p33':prob[32],'p34':prob[33],'p35':prob[34],\
                    'p36':prob[35],'p37':prob[36],'p38':prob[37],'p39':prob[38],'p40':prob[39],\
                    'p41':prob[40],'p42':prob[41],'p43':prob[42],'p44':prob[43],'p45':prob[44],\
                    'p46':prob[45],'p47':prob[46],'p48':prob[47],'p49':prob[48],'p50':prob[49],\
                    'p51':prob[50],'p52':prob[51],'p53':prob[52],'p54':prob[53],'p55':prob[54],\
                    'p56':prob[55],'p57':prob[56],'p58':prob[57],'p59':prob[58],'p60':prob[59],\
                    'p61':prob[60],'p62':prob[61],'p63':prob[62],'p64':prob[63],'dis':dis})
    return out


def MeHperwindow(pat,start,dis,chrom,D,w,ML,depth,optional,MeH=2,dist=1,strand='f'): 
    count=np.zeros((2**w,1))
    m=np.shape(pat)[0]
    #print(count)
    #for i in range(2**w): 
    #    c = 0
    #    for j in range(m):
    #        if (all_pos[i,:]==pat.iloc[j,:]).sum()==w:
    #            c+=1
    #    count[i]=c
    pat=np.array(pat)
    if w==2:
        pat = Counter([str(i[0])+str(i[1]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['00','10','01','11']])
    if w==3:
        pat = Counter([str(i[0])+str(i[1])+str(i[2]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['000','100','010','110','001','101','011','111']])
    if w==4:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['0000','1000','0100','1100','0010','1010','0110','1110','0001',\
                                '1001','0101','1101','0011','1011','0111','1111']])
    if w==5:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[4]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['00000','10000','01000','11000','00100','10100','01100','11100','00010',\
                                '10010','01010','11010','00110','10110','01110','11110','00001','10001','01001','11001','00101',\
                                '10101','01101','11101','00011','10011','01011','11011','00111','10111','01111','11111']])
    if w==6:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[4])+str(i[5]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['000000','100000','010000','110000','001000','101000','011000','111000','000100',\
                                '100100','010100','110100','001100','101100','011100','111100','000010','100010','010010','110010','001010',\
                                '101010','011010','111010','000110', '100110','010110','110110','001110','101110','011110','111110',\
                                '000001','100001','010001','110001','001001','101001','011001','111001','000101',\
                                '100101','010101','110101','001101','101101','011101','111101','000011','100011','010011','110011','001011',\
                                '101011','011011','111011','000111', '100111','010111','110111','001111','101111','011111','111111']])
    
    if MeH==1:  # Abundance based
        score=(((count/m)**2).sum(axis=0))**(-1)
    elif MeH==2: # PWS based
        interaction=np.multiply.outer(count/m,count/m).reshape((2**w,2**w))
        Q=sum(sum(D*interaction))
        #print("Q =",Q)
        if Q==0:
            score=0
        else:
            score=(sum(sum(D*(interaction**2)))/(Q**2))**(-0.5)
    elif MeH==3: #Phylogeny based
        count=count.reshape(2**w)
        count=np.concatenate((count[[0]],count))
        if dist==1 and w==4:
            phylotree=np.append(np.append(np.append(np.append([0],np.repeat(0.5,16)),np.repeat(0.25,6)),[0.5]),np.repeat(0.25,6))
            #phylotree=np.repeat(0,1).append(np.repeat(0.5,16)).append(np.repeat(0.25,6)).append(0.5).append(np.repeat(0.25,6))
            countn=np.zeros(30)
            #count<-rep(0,29)
            countn[1:17]=count[[1,9,5,3,2,13,11,10,7,6,4,15,14,12,8,16]]
            countn[17]=countn[4]+countn[7]
            countn[18]=countn[9]+countn[12]
            countn[19]=countn[1]+countn[2]
            countn[20]=countn[3]+countn[6]
            countn[21]=countn[17]+countn[18]
            countn[22]=countn[19]+countn[20]
            countn[23]=countn[21]+countn[22]
            countn[24]=countn[5]+countn[8]
            countn[25]=countn[10]+countn[13]
            countn[26]=countn[24]+countn[25]
            countn[27]=countn[23]+countn[26]
            countn[28]=countn[11]+countn[14]
            countn[29]=countn[27]+countn[28]
            #Q=sum(sum(phylotree*count))    
        if dist==2 and w==4: 
            phylotree=np.append(np.append(np.append(np.append([0],np.repeat(3,16)),np.repeat(1.5,6)),[3.2,0.8]),np.repeat(2,3),np.repeat(1.5,2))
            #phylotree=c(rep(3,16),rep(1.5,6),3.2,0.8,rep(2,3),1.5,1.5)
            countn=np.zeros(30)
            #print(count)
            countn[1:17]=count[[1,9,5,3,2,13,11,10,7,6,4,15,14,12,8,16]]
            countn[17]=countn[1]+countn[2]
            countn[18]=countn[5]+countn[8]
            countn[19]=countn[3]+countn[6]
            countn[20]=countn[10]+countn[13]
            countn[21]=countn[4]+countn[7]
            countn[22]=countn[11]+countn[14]
            countn[23]=countn[17]+countn[18]
            countn[24]=countn[21]+countn[22]
            countn[25]=countn[19]+countn[20]
            countn[26]=countn[23]+countn[24]
            countn[27]=countn[25]+countn[26]
            countn[28]=countn[9]+countn[12]
            countn[29]=countn[27]+countn[28]
            #Q=sum(phylotree*count)
        if dist==2 and w==3:
            phylotree=np.append(np.append(np.append([0],np.repeat(1.5,8)),np.repeat(0.75,3)),np.repeat(1.5,0.75))
            #phylotree=np.array(0).append(np.repeat(1.5,8)).append(np.repeat(0.75,3)).append(1.5,0.75)
            #phylotree=c(rep(1.5,8),rep(0.75,3),1.5,0.75)
            countn=np.zeros(14)
            countn[1:9]=count[1:9]
            countn[9]=countn[1]+countn[2]
            countn[10]=countn[5]+countn[6]
            countn[11]=countn[3]+countn[4]
            countn[12]=countn[9]+countn[10]
            countn[13]=countn[11]+countn[12]
            #Q=sum(phylotree*count)
        if dist==1 and w==3:
            phylotree=np.append(np.append(np.append([0],np.repeat(0.5,8)),np.repeat(0.25,3)),[0.5,0.25])
            #phylotree=np.array(0).append(np.repeat(0.5,8)).append(np.repeat(0.25,3)).append(0.5,0.25)
            countn=np.zeros(14)
            countn[1:9]=count[1:9]
            countn[9]=countn[1]+countn[2]
            countn[10]=countn[5]+countn[6]
            countn[11]=countn[3]+countn[4]
            countn[12]=countn[9]+countn[10]
            countn[13]=countn[11]+countn[12]
            #print("count = ",count)
            #print("phylotree = ",phylotree)
        Q=sum(phylotree*countn)
        #div=sum(phylotree*((count/Q)**q))**(1/(1-q))
        score=sum(phylotree*((countn/Q)**2))**(-1)
    elif MeH==4: #Entropy
        score=0
        for i in count:
            if i>0:
                score-=(i/m)*np.log2(i/m)/w
    elif MeH==5: #Epipoly
        score=1-((count/m)**2).sum(axis=0)
    out=pd.DataFrame({'chrom':chrom,'pos':start,'MeH':round(div,5),'dis':dis,'ML':round(ML,3),'depth':depth,'strand':strand}, index=[0])    
    
    if optional:
        if d==3:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[0],'p02':count[1],'p03':count[2],'p04':count[3],\
                        'p05':count[4],'p06':count[5],'p07':count[6],'p08':count[7],'MeH':round(div,5),'dis':dis,'ML':round(ML,3),'depth':depth,'strand':strand}, index=[0])     
        if d==4:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[0],'p02':count[1],'p03':count[2],'p04':count[3],\
                        'p05':count[4],'p06':count[5],'p07':count[6],'p08':count[7],'p09':count[8],'p10':count[9],\
                        'p11':count[10],'p12':count[11],'p13':count[12],'p14':count[13],'p15':count[14],\
                        'p16':count[15],'MeH':round(div,5),'dis':dis,'ML':round(ML,3),'depth':depth,'strand':strand}, index=[0])   
        if d==5:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[0],'p02':count[1],'p03':count[2],'p04':count[3],\
                        'p05':count[4],'p06':count[5],'p07':count[6],'p08':count[7],'p09':count[8],'p10':count[9],\
                        'p11':count[10],'p12':count[11],'p13':count[12],'p14':count[13],'p15':count[14],\
                        'p16':count[15],'p17':count[16],'p18':count[17],'p19':count[18],'p20':count[19],\
                        'p21':count[20],'p22':count[21],'p23':count[22],'p24':count[23],'p25':count[24],\
                        'p26':count[25],'p27':count[26],'p28':count[27],'p29':count[28],'p30':count[29],\
                        'p31':count[30],'p32':count[31],'MeH':round(div,5),'dis':dis,'ML':round(ML,3),'depth':depth,'strand':strand}, index=[0])    
        if d==6:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[0],'p02':count[1],'p03':count[2],'p04':count[3],\
                        'p05':count[4],'p06':count[5],'p07':count[6],'p08':count[7],'p09':count[8],'p10':count[9],\
                        'p11':count[10],'p12':count[11],'p13':count[12],'p14':count[13],'p15':count[14],\
                        'p16':count[15],'p17':count[16],'p18':count[17],'p19':count[18],'p20':count[19],\
                        'p21':count[20],'p22':count[21],'p23':count[22],'p24':count[23],'p25':count[24],\
                        'p26':count[25],'p27':count[26],'p28':count[27],'p29':count[28],'p30':count[29],\
                        'p31':count[30],'p32':count[31],'p33':count[32],'p34':count[33],'p35':count[34],\
                        'p36':count[35],'p37':count[36],'p38':count[37],'p39':count[38],'p40':count[39],\
                        'p41':count[40],'p42':count[41],'p43':count[42],'p44':count[43],'p45':count[44],\
                        'p46':count[45],'p47':count[46],'p48':count[47],'p49':count[48],'p50':count[49],\
                        'p51':count[50],'p52':count[51],'p53':count[52],'p54':count[53],'p55':count[54],\
                        'p56':count[55],'p57':count[56],'p58':count[57],'p59':count[58],'p60':count[59],\
                        'p61':count[60],'p62':count[61],'p63':count[62],'p64':count[63],'MeH':round(div,5),'dis':dis,'ML':round(ML,3),'depth':depth,'strand':strand}, index=[0])    
        return out, opt
    else:
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


def CGgenome_scr(bamfile,w,fa,optional,silence=False,dist=1,MeH=2):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
        
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','ML','depth','strand'])
    
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','ML','depth','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','ML','depth','strand'])    

    
    neverr = never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) #1:Hamming distance
    
    start=datetime.datetime.now()
    
    for pileupcolumn in samfile.pileup():
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        # Forward
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
                
        # Reverse
        if pileupcolumn.pos>1:
            if (fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos+1)=='GC'):        
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                        dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        dfr2 = pd.DataFrame(data=dr)
                        #df2.head()
                        tempr=tempr.append(dfr2, ignore_index=True)
                #temp.head()
                if (not tempr.empty):
                    #temp.head()
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()

        # Impute and estimate
        if never and aggreC.shape[1] == (2*w):
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='f')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'f'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
            
        # Reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='r')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='f')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'f'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CG. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
        
        # reverse
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='r')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
                            if optional:
                                Resultopt.to_csv(r"MeHdata/CG_opt_%s.csv"%(filename),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CG. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CG_opt_%s.csv"%(filename),index = False, header=True)
                            
    print("Done. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))
            
    #samfile.close()  
    
def CHHgenome_scr(bamfile,w,fa,optional,silence=False,dist=1,MeH=2):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
        
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','ML','depth','strand'])
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','ML','depth','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','ML','depth','strand'])  
    neverr = never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) #1:Hamming distance
    
    start=datetime.datetime.now()
    
    for pileupcolumn in samfile.pileup():
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CHH %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        # forward
        if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)!='G':        
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
                
        # reverse
        if pileupcolumn.pos>2:
            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)!='C':        
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                        d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2 = pd.DataFrame(data=d)
                        #df2.head()
                        tempr=tempr.append(df2, ignore_index=True)
                #temp.head()
                if (not tempr.empty):
                    #temp.head()
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()   

        if never and aggreC.shape[1] == (2*w):
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    # MeH eligibility
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='f')
                        ResultPW=ResultPW.append(toappend)
                    # else output methylation level
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'f'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
            
        # reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                    
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    #if enough_reads(window,w,complete=True):

                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='r')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)
                    
            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    # MeH eligibility
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='f',optional=optional)
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='f',optional=optional)
                        ResultPW=ResultPW.append(toappend)
                        
    
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'f'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
                            if optional:
                                Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
                            
                            if not silence: 
                                print("Checkpoint CHH. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
            
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='r',optional=optional)
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,strand='r',optional=optional)
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
                            if optional:
                                Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHH. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w 
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
                            
    print("Done CHH. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))
            
def CHGgenome_scr(bamfile,w,fa,optional,silence=False,dist=1,MeH=2):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
        
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','ML','depth','strand'])
    
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','ML','depth','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','ML','depth','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','ML','depth','strand'])    

    neverr = never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) #1:Hamming distance
    
    start=datetime.datetime.now()
    
    for pileupcolumn in samfile.pileup():
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CHG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)!='G':        
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
        
        # reverse
        if pileupcolumn.pos>2:
            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)=='C':        
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # G
                        dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2r = pd.DataFrame(data=dr)
                        #df2.head()
                        tempr=tempr.append(df2r, ignore_index=True)
                #temp.head()
                if (not tempr.empty):
                    #temp.head()
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()        

        if never and aggreC.shape[1] == (2*w):
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    #if enough_reads(window,w,complete=True):                 
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='f')
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='f')
                        ResultPW=ResultPW.append(toappend)
                        
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'f'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

                    
            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
        
        # reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i<w:
                    window = meth.iloc[:,range(i,i+w)].values
                    #if enough_reads(window,w,complete=True):
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='r')
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='r')
                        ResultPW=ResultPW.append(toappend)
                        
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            #MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                            #                dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                            #                chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth)
                            ResultPW=ResultPW.append(toappend)

            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
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
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                # if eligible for imputation
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='f')
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='f')
                        ResultPW=ResultPW.append(toappend)
                        
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth}, index=[0])    
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
                            if optional:
                                Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHG. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
        # reverse
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['C','G'],1)
            aggreR = aggreR.replace(['A','T'],0)
            aggreR = aggreR.replace(['N'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            methtemp = meth.copy()
            # impute once if valid
            for i in range(0,meth.shape[1]-w+1,1):
                window = meth.iloc[:,range(i,i+w)].values
                # if eligible for imputation
                MC=(window==1).sum(axis=0)[0]
                UC=(window==0).sum(axis=0)[0]
                depth=MC+UC
                if enough_reads(window,w,complete=False):
                    window=pd.DataFrame(data=impute(window,w))
                    ind=np.where(window.notnull().sum(axis=1)==w)[0]
                    methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
            meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(0,meth.shape[1]-w+1,1):
                if i>w-2 and i<2*w:
                    window = meth.iloc[:,range(i,i+w)].values
                    #if enough_reads(window,w,complete=True):
                    if enough_reads(window,w,complete=True):
                        ML=float(MC)/float(depth)
                        matforMH=getcomplete(window,w)
                        if optional:
                            toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='r')
                            Resultopt=Resultopt.append(opt)
                        else:
                            toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=2,ML=ML,depth=depth,optional=optional,strand='r')
                        ResultPW=ResultPW.append(toappend)
                    else:
                        if depth>3:
                            ML=float(MC)/float(depth)
                            toappend=pd.DataFrame({'chrom':chrom,'pos':meth.iloc[:,range(i,i+w)].columns[0],'MeH':np.nan,'dis':np.nan,'ML':round(ML,3),'depth':depth,'strand':'r'}, index=[0])    
                            ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
                            if optional:
                                Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHG. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True) 
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                            
    
    print("Done CHG. For sample %s %s: %s results obtained up to position %s." % (filename,chrom,ResultPW.shape[0],pileupcolumn.pos))
            
        
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
        if(outfile_Size >=337374182 and sum_Outfile_Size <= bamsize):
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
parser.add_argument("-w", "--windowsize",type=int, default=4 ,help='number of CGs')
#parser.add_argument("-c", "--chromosome",type=str,help='chromosome')
parser.add_argument("-c", "--cores",type=int, default=4, help='number of cores')
parser.add_argument("-m", "--MeH",type=int, default=2, help='Methylation heterogeneity score 1:Abundance 2:PW 3:Phylogeny')
parser.add_argument("-d", "--dist",type=int, default=1, help='Distance between methylation patterns 1:Hamming 2:WDK')
#parser.add_argument("-g", "--context",type=int, default=1, help='Cytosine context 1:CG 2:CHG 3:CHH')
parser.add_argument('--CG', default=False, action='store_true')
parser.add_argument('--CHG', default=False, action='store_true')
parser.add_argument('--CHH', default=False, action='store_true')
parser.add_argument('--opt', default=False, action='store_true')



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
    
    #if 'cores' in args: 
    #    num_cores = args.cores
    #else:
    #    num_cores = 4
        
    Parallel(n_jobs=args.cores)(delayed(split_bam)(bamfile,Folder=Folder) for bamfile in bam_list)
    
    spbam_list = []
    tempfiles = os.listdir(Folder)
    for file in tempfiles:
        filename, file_extension = os.path.splitext(file)
        if file_extension=='.bam' and filename not in bam_list:
            spbam_list.append(filename)
    #print(spbam_list)
        
    #start=t.time()
    if args.CG:
        con='CG'
        Parallel(n_jobs=args.cores)(delayed(CGgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt) for bamfile in spbam_list)
    
        # merge .csv within sample
        for file in spbam_list:
            #filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' +file + '.csv'
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
        if args.opt:
            for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_opt_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_opt_' +file + '.csv'
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
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                Tomerge = Tomerge.drop(columns=['dis','ML','depth'])
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','pos','strand'])
                Result = Result.drop_duplicates() 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result = Result.drop(columns=['dis','ML','depth'])
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)

        Result.to_csv(Folder + con + '_' +'Results.csv' ,index = False,header=True)
        print("All done. ",len(bam_list)," bam files processed and merged for CG.")
    
    if args.CHG:
        con='CHG'
        Parallel(n_jobs=args.cores)(delayed(CHGgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt) for bamfile in spbam_list)
    
        # merge .csv within sample
        for file in spbam_list:
            #filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' +file + '.csv'
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
        if args.opt:
            for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_opt_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_opt_' +file + '.csv'
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
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                Tomerge = Tomerge.drop(columns=['dis','ML','depth'])
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','pos','strand'])
                Result = Result.drop_duplicates() 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result = Result.drop(columns=['dis','ML','depth'])
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)

        Result.to_csv(Folder + con + '_' +'Results.csv' ,index = False,header=True)
        print("All done. ",len(bam_list)," bam files processed and merged for CHG.")

    if args.CHH:
        con='CHH'
        Parallel(n_jobs=args.cores)(delayed(CHHgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt) for bamfile in spbam_list)
    
        # merge .csv within sample
        for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' +file + '.csv'
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
        
        if args.opt:
            for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("sample = ",sample)
            if not sample == filename:
                res_dir = Folder + con + '_opt_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_opt_' +file + '.csv'
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
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                Tomerge = Tomerge.drop(columns=['dis','ML','depth'])
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','pos','strand'])
                Result = Result.drop_duplicates() 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result = Result.drop(columns=['dis','ML','depth'])
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)

        Result.to_csv(Folder + con + '_' +'Results.csv' ,index = False,header=True)

    for filename in spbam_list:
        file = Folder + filename + '.bam'
        os.remove(file)
        
    print("All done. ",len(bam_list)," bam files processed and merged for CHH.")

# FINAL
