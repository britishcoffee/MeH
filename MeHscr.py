# Mar21, 2022

##
#---------------------------------------------------------------------
# SERVER only input all files (.bam and .fa) output MeH matrix in .csv
# August 3, 2021 clean
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

    
def open_log(fname):
    open_log.logfile = open(fname, 'w', 1)
    

def logm(message):
    log_message = "[%s] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print(log_message),
    open_log.logfile.write(log_message)

def close_log():
    open_log.logfile.close()
    
    
# Count # of windows with enough reads for complete/impute
def coverage(methbin, complete, w):
    meth = methbin.iloc[:, methbin.columns != 'Qname']
    if len(meth.columns) >= w:
        tot = len(meth.columns) - w + 1
        count = sum(enough_reads(window=meth.iloc[:, i:i+w].copy(), complete=complete, w=w) for i in range(tot))
        return count / tot * 100
    else:
        return 0
    

# Check whether a window has enough reads for complete/impute
def enough_reads(window,w,complete):
    temp=np.isnan(window).sum(axis=1)==0
    if complete: # For heterogeneity estimation
        return temp.sum()>=2**w
    else:  # for imputation
        tempw1=np.isnan(window).sum(axis=1)==1
        return temp.sum()>=2**(w-2) and tempw1.sum()>0
    

def impute(window, w):
    full_ind = np.where(np.isnan(window).sum(axis=1) == 0)[0]
    part_ind = np.where(np.isnan(window).sum(axis=1) == 1)[0]

    for i in range(len(part_ind)):
        pos = np.where(np.isnan(window[part_ind[i], :]))[0]
        if np.unique(window[np.where(~np.isnan(window[:, pos]))[0], pos]).shape[0] == 1:
            window[part_ind[i], pos] = window[np.where(~np.isnan(window[:, pos]))[0], pos][0]
        else:
            sam = np.where(np.sum(window[part_ind[i], :] == window[full_ind, :], axis=1) == w - 1)[0]
            if len(sam) > 0:
                s1 = random.choice(sam)
                s = window[full_ind[s1], pos]
            else:
                valid_values = window[np.where(~np.isnan(window[:, pos]))[0], pos]
                s = random.choice(valid_values.tolist())
            window[part_ind[i], pos] = np.float64(s)
    
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
def window_summ(pat, start, dis, chrom):
    m, d = pat.shape
    all_pos = np.unpackbits(np.arange(2**d, dtype=np.uint8)[:, np.newaxis], axis=1)[:, -d:]
    
    counts = (np.logical_and(np.all(all_pos[:, np.newaxis, :] == pat.values, axis=2), np.sum(~pat.isnull().values, axis=1) == d)).sum(axis=1)
    
    column_names = ['chrom', 'pos'] + ['p{:02d}'.format(i + 1) for i in range(2 ** d)] + ['dis']
    out_dict = {'chrom': chrom, 'pos': start, 'dis': dis}
    out_dict.update(zip(column_names[2:-1], counts))
    
    out = pd.DataFrame(out_dict)
    return out


def MeHperwindow(pat,start,dis,chrom,D,w,optional,MeH=2,dist=1,strand='f'): 
    count=np.zeros((2**w,1))
    m=np.shape(pat)[0]
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
        score=sum(phylotree*((countn/Q)**2))**(-1)
    elif MeH==4: #Entropy
        score=0
        for i in count:
            if i>0:
                score-=(i/m)*np.log2(i/m)/w
    elif MeH==5: #Epipoly
        score=1-((count/m)**2).sum(axis=0)
    
    if optional:
        if MeH!=3:
            count=count.reshape(2**w)
            count=np.concatenate((count[[0]],count))
        if w==3:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'MeH':round(score,5),'dis':dis,'strand':strand}, index=[0])     
        if w==4:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'MeH':round(score,5),'dis':dis,'strand':strand}, index=[0])   
        if w==5:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'p17':count[17],'p18':count[18],'p19':count[19],'p20':count[20],\
                        'p21':count[21],'p22':count[22],'p23':count[23],'p24':count[24],'p25':count[25],\
                        'p26':count[26],'p27':count[27],'p28':count[28],'p29':count[29],'p30':count[30],\
                        'p31':count[31],'p32':count[32],'MeH':round(score,5),'dis':dis,'strand':strand}, index=[0])    
        if w==6:
            opt=pd.DataFrame({'chrom':chrom,'pos':start,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'p17':count[17],'p18':count[18],'p19':count[19],'p20':count[20],\
                        'p21':count[21],'p22':count[22],'p23':count[23],'p24':count[24],'p25':count[25],\
                        'p26':count[26],'p27':count[27],'p28':count[28],'p29':count[29],'p30':count[30],\
                        'p31':count[31],'p32':count[32],'p33':count[33],'p34':count[34],'p35':count[35],\
                        'p36':count[36],'p37':count[37],'p38':count[38],'p39':count[39],'p40':count[40],\
                        'p41':count[41],'p42':count[42],'p43':count[43],'p44':count[44],'p45':count[45],\
                        'p46':count[46],'p47':count[47],'p48':count[48],'p49':count[49],'p50':count[50],\
                        'p51':count[51],'p52':count[52],'p53':count[53],'p54':count[54],'p55':count[55],\
                        'p56':count[56],'p57':count[57],'p58':count[58],'p59':count[59],'p60':count[60],\
                        'p61':count[61],'p62':count[62],'p63':count[63],'p64':count[64],'MeH':round(score,5),'dis':dis,'strand':strand}, index=[0])    
        return out, opt
    else:
        out=pd.DataFrame({'chrom':chrom,'pos':start,'MeH':round(score,5),'dis':dis,'strand':strand}, index=[0])    
        return out


def CGgenome_scr(bamfile,w,fa,optional,melv,silence=False,dist=1,MeH=2,imp=True):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    coverage = cov_context = 0
    # load bamfile
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    # load reference genome
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    
    # initialise data frame for genome screening (load C from bam file)
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    
    # initialise data frame for output 
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','strand'])
    if melv:
        ResML = pd.DataFrame(columns=['chrom','pos','ML','strand','depth'])
    
    
    # if user wants to output compositions of methylation patterns at every eligible window, initialise data frame
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','strand'])
        if w==6:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','strand'])    

    
    neverr = never = True
    
    # all methylation patterns for Methylation heterogeneity evaluation
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    # distance matrix, also for Methylation heterogeneity evaluation
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) # 1:Hamming distance, 2: WDK
    
    start=datetime.datetime.now()
    
    # vector for saving methylation statuses before imputation
    MU=np.zeros((2,w))
    
    # screen bamfile by column
    for pileupcolumn in samfile.pileup():
        coverage += 1 
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        # Forward strand, check if 'CG' in reference genome 
        if (fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+2)=='CG'):        
            cov_context += 1
            temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
            pileupcolumn.set_min_base_quality(0)
            # append reads in the column
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                    d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                    df2 = pd.DataFrame(data=d)
                    temp=temp.append(df2, ignore_index=True)
            if melv:
                temp2 = temp.replace(['C'],1)
                temp2 = temp2.replace(['G'],0)
                temp2 = temp2.replace(['A','T','N'],np.nan)
                temp2 = temp2.drop('Qname',axis=1)
                MC=(temp2==1).sum(axis=0).to_numpy()
                UC=(temp2==0).sum(axis=0).to_numpy()
                depth=MC+UC
                if depth>3:
                    toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                        'strand':'f','depth':depth,'ML':float(MC)/depth}, index=[0])
                    ResML=ResML.append(toappend)
            # merge with other columns
            if (not temp.empty):
                aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                aggreC = aggreC.drop_duplicates()
                
        # Reverse strand, check if 'CG' in reference genome 
        if pileupcolumn.pos>1:
            if (fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos+1)=='CG'):
                cov_context += 1
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                        dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        dfr2 = pd.DataFrame(data=dr)
                        tempr=tempr.append(dfr2, ignore_index=True)
                if melv:
                    temp2 = tempr.replace(['G'],1)
                    temp2 = temp2.replace(['C'],0)
                    temp2 = temp2.replace(['A','T','N'],np.nan)
                    temp2 = temp2.drop('Qname',axis=1)
                    MC=(temp2==1).sum(axis=0).to_numpy()
                    UC=(temp2==0).sum(axis=0).to_numpy()
                    depth=MC+UC
                    if depth>3:
                        toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                        'strand':'r','depth':depth,'ML':float(MC)/depth}, index=[0])
                        ResML=ResML.append(toappend)
                        
                if (not tempr.empty):
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()

        # Impute and estimate, if there are 2w-1 columns
        if never and aggreC.shape[1] == (2*w):
            # C/G to 1, rest to 0, N to NA
            never = False
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['A','N','G'],np.nan)
            methbin = aggreC 
            meth = methbin.copy()
            # remove read ID
            meth = meth.drop('Qname',axis=1)
            # back up for imputation
            if imp:
                methtemp = meth.copy()
                # imputation by sliding window of 1 C
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # save methylation statuses before imputation
                    # check if eligible for imputation, impute
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                # overwrite imputed window
                meth = methtemp.copy()
            # Evaluate methylation level and methylation heterogeneity and append to result
            for i in range(0,w,1): # w windows
                window = meth.iloc[:,range(i,i+w)].values
                # check if enough complete patterns for evaluating MeH
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    # if need to output methylation patterns
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                        Resultopt=Resultopt.append(opt)
                        # evaluate and output MeH 
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                    ResultPW=ResultPW.append(toappend)
                    
            # remove 1 column
            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            # drop rows with no values
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
            
        # Reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['C','N','T'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(0,w,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                    ResultPW=ResultPW.append(toappend)
            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
        #------------------
        #  SECONDARY CASE
        #------------------

        if (aggreC.shape[1] == (3*w-1)):
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['A','N','G'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
            # compute coverage and output summary
            for i in range(w-1,2*w-1,1):
            #for i in range(0,meth.shape[1]-w+1,1):
                #if i>w-2 and i<2*w:
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                    ResultPW=ResultPW.append(toappend)

                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CG_ML_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CG_opt_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))
            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
        
        # reverse
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['C','N','T'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(w-1,2*w-1,1):
                window = meth.iloc[:,range(i,i+w)].values   
                if enough_reads(window,w,complete=True): 
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                    ResultPW=ResultPW.append(toappend)
                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CG_ML_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CG_opt_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CG_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CG_opt_%s.csv"%(filename),index = False, header=True)
        if melv:
            ResML.to_csv(r"MeHdata/CG_ML_%s.csv"%(filename),index = False, header=True)
                               
    return sample, coverage, cov_context, 'CG'        
    print("Done CG for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))
            
    #samfile.close()  
    
def CHHgenome_scr(bamfile,w,fa,optional,melv,silence=False,dist=1,MeH=2,imp=True):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    coverage = cov_context = 0
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','strand'])
    if melv:
        ResML = pd.DataFrame(columns=['chrom','pos','ML','depth','strand'])

                            
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','strand'])  
    neverr = never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) #1:Hamming distance
    
    start=datetime.datetime.now()
    MU=np.zeros((2,w))
    
    for pileupcolumn in samfile.pileup():
        coverage += 1
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CHH %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        # forward
        if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)!='G':        
            cov_context += 1
            temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
            pileupcolumn.set_min_base_quality(0)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                    d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                    df2 = pd.DataFrame(data=d)
                    #df2.head()
                    temp=temp.append(df2, ignore_index=True)
            #temp.head()
            if melv:
                temp2 = temp.replace(['C'],1)
                temp2 = temp2.replace(['T'],0)
                temp2 = temp2.replace(['A','G','N'],np.nan)
                temp2 = temp2.drop('Qname',axis=1)
                MC=(temp2==1).sum(axis=0).to_numpy()
                UC=(temp2==0).sum(axis=0).to_numpy()
                depth=MC+UC
                if depth>3:
                    toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                        'strand':'f','depth':depth,'ML':float(MC)/depth}, index=[0])
                    ResML=ResML.append(toappend)
            if (not temp.empty):
                #temp.head()
                aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                aggreC = aggreC.drop_duplicates()
                
        # reverse
        if pileupcolumn.pos>2:
            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)!='C':        
                cov_context += 1
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                        d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2 = pd.DataFrame(data=d)
                        #df2.head()
                        tempr=tempr.append(df2, ignore_index=True)
                #temp.head()
                if melv:
                    temp2 = tempr.replace(['G'],1)
                    temp2 = temp2.replace(['A'],0)
                    temp2 = temp2.replace(['C','T','N'],np.nan)
                    temp2 = temp2.drop('Qname',axis=1)
                    MC=(temp2==1).sum(axis=0).to_numpy()
                    UC=(temp2==0).sum(axis=0).to_numpy()
                    depth=MC+UC
                    if depth>3:
                        toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                        'strand':'r','depth':depth,'ML':float(MC)/depth}, index=[0])
                        ResML=ResML.append(toappend)
                if (not tempr.empty):
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()   

        if never and aggreC.shape[1] == (2*w):
            never = False
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['A','N','G'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(0,w,1):
                window = meth.iloc[:,range(i,i+w)].values
                # MeH eligibility
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                    ResultPW=ResultPW.append(toappend)
                        
            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
            
        # reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['C','N','T'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(0,w,1):
                window = meth.iloc[:,range(i,i+w)].values
                #if enough_reads(window,w,complete=True):
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                    ResultPW=ResultPW.append(toappend)
            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            
        #------------------
        #  SECONDARY CASE
        #------------------

        if (aggreC.shape[1] == (3*w-1)):
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['N','G','A'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(w-1,2*w-1,1):
                window = meth.iloc[:,range(i,i+w)].values
                # MeH eligibility
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='f',optional=optional)
                    ResultPW=ResultPW.append(toappend)
                        
                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CHH_ML_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CHH. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
            
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['N','T','C'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(w-1,2*w-1,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,strand='r',optional=optional)
                    ResultPW=ResultPW.append(toappend)
                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CHH_ML_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CHH. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w 
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CHH_%s.csv"%(filename),index = False, header=True)
        if melv:
            ResML.to_csv(r"MeHdata/CHH_ML_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CHH_opt_%s.csv"%(filename),index = False, header=True)
    return sample, coverage, cov_context, 'CHH'                        
    print("Done CHH for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))
            
def CHGgenome_scr(bamfile,w,fa,optional,melv,silence=False,dist=1,MeH=2,imp=True):
    filename, file_extension = os.path.splitext(bamfile)
    sample = str.split(filename,'_')[0]
    #directory = "Outputs/" + str(sample) + '.csv' #original filename of .bams
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    coverage = cov_context = 0
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    ResultPW = pd.DataFrame(columns=['chrom','pos','MeH','dis','strand'])
    if melv:
        ResML = pd.DataFrame(columns=['chrom','pos','ML','depth','strand'])
        
    if optional:
        if w==3:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08',\
                         'MeH','dis','strand'])
        if w==4:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','MeH','dis','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','MeH','dis','strand'])
        if w==5:
            Resultopt = pd.DataFrame(columns=\
                        ['chrom','pos','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'MeH','dis','strand'])    

    neverr = never = True
    #chr_lengths = fastafile.get_reference_length(chrom)
    all_pos=np.zeros((2**w,w))
    for i in range(w): 
        all_pos[:,i]=np.linspace(0,2**w-1,2**w)%(2**(i+1))//(2**i)
    
    D=PattoDis(pd.DataFrame(all_pos),dist=dist) #1:Hamming distance
    MU=np.zeros((2,w))
    start=datetime.datetime.now()
    
    for pileupcolumn in samfile.pileup():
        coverage += 1
        chrom = pileupcolumn.reference_name
        if not silence:
            if (pileupcolumn.pos % 2000000 == 1):
                print("CHG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))
        
        if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)=='G':        
            cov_context += 1
            temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
            pileupcolumn.set_min_base_quality(0)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                    d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                    df2 = pd.DataFrame(data=d)
                    #df2.head()
                    temp=temp.append(df2, ignore_index=True)
            #temp.head()
            if melv:
                temp2 = temp.replace(['C'],1)
                temp2 = temp2.replace(['T'],0)
                temp2 = temp2.replace(['A','G'],np.nan)
                temp2 = temp2.drop('Qname',axis=1)
                MC=(temp2==1).sum(axis=0).to_numpy()
                UC=(temp2==0).sum(axis=0).to_numpy()
                depth=MC+UC
                if depth>3:
                    toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                            'strand':'f','depth':depth,'ML':float(MC)/float(MC+UC)}, index=[0])
                    ResML=ResML.append(toappend)
            if (not temp.empty):
                #temp.head()
                aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                aggreC = aggreC.drop_duplicates()
        
        # reverse
        if pileupcolumn.pos>2:
            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)=='C':        
                cov_context += 1
                tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # G
                        dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2r = pd.DataFrame(data=dr)
                        #df2.head()
                        tempr=tempr.append(df2r, ignore_index=True)
                #temp.head()
                if melv:
                    temp2 = tempr.replace(['G'],1)
                    temp2 = temp2.replace(['A'],0)
                    temp2 = temp2.replace(['C','T'],np.nan)
                    temp2 = temp2.drop('Qname',axis=1)
                    MC=(temp2==1).sum(axis=0).to_numpy()
                    UC=(temp2==0).sum(axis=0).to_numpy()
                    depth=MC+UC
                    if depth>3:
                        toappend=pd.DataFrame({'chrom':chrom,'pos':temp2.columns[0], \
                                            'strand':'r','depth':depth,'ML':float(MC)/float(MC+UC)}, index=[0])
                        ResML=ResML.append(toappend)
                if (not tempr.empty):
                    #temp.head()
                    aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                    aggreR = aggreR.drop_duplicates()        

        if never and aggreC.shape[1] == (2*w):
            never = False
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['A','G','N'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(0,w,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='f')
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='f')
                    ResultPW=ResultPW.append(toappend)
               
            aggreC = aggreC.drop(meth.columns[0:1],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
        
        # reverse
        if neverr and aggreR.shape[1] == (2*w):
            neverr = False
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['N','C','T'],np.nan)
            methbin = aggreR # backup
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(0,w,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='r')
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='r')
                    ResultPW=ResultPW.append(toappend)

            aggreR = aggreR.drop(meth.columns[0:1],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True)
            #total += w
        #------------------
        #  SECONDARY CASE
        #------------------

        if (aggreC.shape[1] == (3*w-1)):
            aggreC = aggreC.replace(['C'],1)
            aggreC = aggreC.replace(['T'],0)
            aggreC = aggreC.replace(['N','A','G'],np.nan)
            methbin = aggreC # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(w-1,2*w-1,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='f')
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                    dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                    chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='f')
                    ResultPW=ResultPW.append(toappend)

                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CHG_ML_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CHG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

            aggreC = aggreC.drop(meth.columns[0:w],axis=1)
            aggreC.dropna(axis = 0, thresh=2, inplace = True)
            #print(aggreC)
            #total += w
        # reverse
        if (aggreR.shape[1] == (3*w-1)):
            aggreR = aggreR.replace(['G'],1)
            aggreR = aggreR.replace(['A'],0)
            aggreR = aggreR.replace(['N','T','C'],np.nan)
            methbin = aggreR # backup
            #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
            meth = methbin.copy()
            meth = meth.drop('Qname',axis=1)
            if imp:
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
            for i in range(w-1,2*w-1,1):
                window = meth.iloc[:,range(i,i+w)].values
                if enough_reads(window,w,complete=True):
                    matforMH=getcomplete(window,w)
                    if optional:
                        toappend,opt=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='r')
                        Resultopt=Resultopt.append(opt)
                    else:
                        toappend=MeHperwindow(pd.DataFrame(matforMH),start=meth.iloc[:,range(i,i+w)].columns[0],\
                                        dis=meth.iloc[:,range(i,i+w)].columns[w-1]-meth.iloc[:,range(i,i+w)].columns[0],\
                                        chrom=chrom,D=D,w=w,dist=dist,MeH=MeH,optional=optional,strand='r')
                    ResultPW=ResultPW.append(toappend)

                    if ResultPW.shape[0] % 100000 == 1:   
                        ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
                        if optional:
                            Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                        if melv:
                            ResML.to_csv(r"MeHdata/CHG_ML_%s.csv"%(filename),index = False, header=True)
                        if not silence: 
                            print("Checkpoint CHG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos+1))

            aggreR = aggreR.drop(meth.columns[0:w],axis=1)
            aggreR.dropna(axis = 0, thresh=2, inplace = True) 
            
    if ResultPW.shape[0]>0:   
        ResultPW.to_csv(r"MeHdata/CHG_%s.csv"%(filename),index = False, header=True)
        if melv:
            ResML.to_csv(r"MeHdata/CHG_ML_%s.csv"%(filename),index = False, header=True)
        if optional:
            Resultopt.to_csv(r"MeHdata/CHG_opt_%s.csv"%(filename),index = False, header=True)
                            
    return sample, coverage, cov_context, 'CHG'
    print("Done CHG for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos+1))
            
        
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
    print("fileout",fileout)
    header = samfile.header
        
    outfile = pysam.Samfile(fileout, "wb", header = header)
    sum_Outfile_Size=0
    for reads in samfile.fetch():
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
            print("fileout",fileout)
            outfile = pysam.Samfile(fileout, "wb",header = header)
            
    outfile.close()
    pysam.index(fileout)


    

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-w", "--windowsize",type=int, default=4 ,help='number of CGs')
parser.add_argument("-c", "--cores",type=int, default=4, help='number of cores')
parser.add_argument("-m", "--MeH",type=int, default=2, help='Methylation heterogeneity score 1:Abundance 2:PW 3:Phylogeny')
parser.add_argument("-d", "--dist",type=int, default=1, help='Distance between methylation patterns 1:Hamming 2:WDK')
parser.add_argument("--CG", default=False, action='store_true', help='Include genomic context CG')
parser.add_argument("--CHG", default=False, action='store_true', help='Include genomic context CHG')
parser.add_argument("--CHH", default=False, action='store_true', help='Include genomic context CHH')
parser.add_argument("--opt", default=False, action='store_true', help='Outputs compositions of methylation patterns')
parser.add_argument('--mlv', default=False, action='store_true', help='Outputs methylation levels')
parser.add_argument('--imp', default=False, action='store_true', help='Whether to implement BSImp (impute if valid)')



args = parser.parse_args()

import sys
import time
import os
import pandas as pd
import multiprocessing
from joblib import Parallel, delayed

#num_cores = multiprocessing.cpu_count()
                                                
if __name__ == "__main__":
    
    open_log('MeHscreening.log')
    logm("Call genome screening.")
    
    #start = time.time()
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
    
    topp = pd.DataFrame(columns=['sample','coverage','context_coverage','context'])    
    #CG = []
    #start=t.time()
    if args.CG:
        con='CG'
        CG=Parallel(n_jobs=args.cores)(delayed(CGgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt,melv=args.mlv,imp=args.imp) for bamfile in spbam_list)
        
        logm("Merging MeH within samples for CG.")
        # merge MeH within sample
        for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("Merging within sample",sample,"...")
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    #os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    
        # not into bins of 400bp
        if args.opt:
            for file in spbam_list:
                filename, file_extension = os.path.splitext(file)
                sample = str.split(file,'_')[0]
                #print("sample = ",sample)
                if not sample == filename:
                    res_dir = Folder + con + '_opt_' + str(sample) + '.csv'
                    toapp_dir = Folder + con + '_opt_' +file + '.csv'
                    if os.path.exists(res_dir):
                        Tomod = pd.read_csv(res_dir) 
                        Toappend = pd.read_csv(toapp_dir)
                        Tomod = Tomod.append(Toappend)
                        Tomod.to_csv(res_dir,index = False, header = True)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)

        #os.chdir('../')
        #os.chdir(outputFolder)
        
        logm("Merging ML within samples for CG.")
        # append ML within samples
        if args.mlv:
            for file in spbam_list:
                filename, file_extension = os.path.splitext(file)
                sample = str.split(file,'_')[0]
                res_dir = Folder + con + '_ML_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_ML_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    #os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    #print(Toappend)
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    
        logm("Merging ML between samples for CG.")
        # merge ML between samples
        if args.mlv:
            for sample in bam_list: 
                tomerge_dir = Folder +  con + '_ML_' + str(sample) + '.csv' 
                res_dir = Folder +  con + '_ML_' + 'Results.csv'
                if os.path.exists(res_dir):
                    Result = pd.read_csv(res_dir) 
                    Tomerge = pd.read_csv(tomerge_dir)
                    Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                    Tomerge = Tomerge.rename(columns={'ML': sample})
                    Result=Result.merge(Tomerge, on=['chrom','bin','strand'])
                    Result.dropna(axis = 0, thresh=4, inplace = True) 
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
                else:
                    Result = pd.read_csv(tomerge_dir)
                    Result = Result.rename(columns={'ML': sample})
                    #Result = Result.drop(columns=['counts','pos','depth','dis'])
                    Result.dropna(axis = 0, thresh=4, inplace = True)
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
        logm("Merging MeH between samples for CG.")              
        # merge MeH between samples
        for sample in bam_list: 
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                
                Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','bin','strand'])
                Result.dropna(axis = 0, thresh=4, inplace = True) 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result.head()
                
                Result.dropna(axis = 0, thresh=4, inplace = True)
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)


        Result.to_csv(Folder + con + '_' +'Results.csv' ,index = False,header=True)
        print("All done.",len(bam_list),"bam files processed and merged for CG.")
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CG.")
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)
        
          
    if args.CHG:
        con='CHG'
        CG=Parallel(n_jobs=args.cores)(delayed(CHGgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt,melv=args.mlv,imp=args.imp) for bamfile in spbam_list)
        
        logm("Merging MeH within samples for CHG.")  
        for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("Merging within sample",sample,"...")
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    
        # not into bins of 400bp
        if args.opt:
            for file in spbam_list:
                filename, file_extension = os.path.splitext(file)
                sample = str.split(file,'_')[0]
                #print("sample = ",sample)
                if not sample == filename:
                    res_dir = Folder + con + '_opt_' + str(sample) + '.csv'
                    toapp_dir = Folder + con + '_opt_' +file + '.csv'
                    if os.path.exists(res_dir):
                        Tomod = pd.read_csv(res_dir) 
                        Toappend = pd.read_csv(toapp_dir)
                        Tomod = Tomod.append(Toappend)
                        Tomod.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)

        logm("Merging ML within samples for CHG.")    
        # append ML within samples
        if args.mlv:
            for file in spbam_list:
                filename, file_extension = os.path.splitext(file)
                sample = str.split(file,'_')[0]
                res_dir = Folder + con + '_ML_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_ML_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    #Count=Count.drop_duplicates()
                    #print(Count)
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    #print(Toappend)
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    os.remove(toapp_dir)
        
        logm("Merging MeH between samples for CHG.")            
        # merge MeH between samples
        for sample in bam_list: 
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                #Tomerge = Tomerge.drop(columns=['dis','ML','depth'])
                Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','bin','strand'])
                Result.dropna(axis = 0, thresh=4, inplace = True) 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result.head()
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
                
        logm("Merging ML between samples for CHG.")     
        # merge ML between samples
        if args.mlv:
            for sample in bam_list: 
                tomerge_dir = Folder +  con + '_ML_' + str(sample) + '.csv' 
                res_dir = Folder +  con + '_ML_' + 'Results.csv'
                if os.path.exists(res_dir):
                    Result = pd.read_csv(res_dir) 
                    Tomerge = pd.read_csv(tomerge_dir)
                    Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                    Tomerge = Tomerge.rename(columns={'ML': sample})
                    Result=Result.merge(Tomerge, on=['chrom','bin','strand'])
                    Result.dropna(axis = 0, thresh=4, inplace = True) 
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
                else:
                    Result = pd.read_csv(tomerge_dir)
                    Result = Result.rename(columns={'ML': sample})
                    #Result = Result.drop(columns=['counts','pos','depth','dis'])
                    Result.dropna(axis = 0, thresh=4, inplace = True)
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
                    
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CHG.")
        
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)
        
    if args.CHH:
        con='CHH'
        CG=Parallel(n_jobs=args.cores)(delayed(CHHgenome_scr)(bamfile,w=args.windowsize,fa=fa,MeH=args.MeH,dist=args.dist,optional=args.opt,melv=args.mlv,imp=args.imp) for bamfile in spbam_list)
    
        logm("Merging MeH within samples for CHH.")
        # merge MeH within sample
        for file in spbam_list:
            filename, file_extension = os.path.splitext(file)
            sample = str.split(file,'_')[0]
            print("Merging within sample",sample,"...")
            if not sample == filename:
                res_dir = Folder + con + '_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    #os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Toappend=Toappend.drop(columns=['pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'MeH': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    
        # not into bins of 400bp
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
                        os.remove(toapp_dir)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)

        logm("Merging ML within samples for CHH.")
        # append ML within samples
        if args.mlv:
            for file in spbam_list:
                filename, file_extension = os.path.splitext(file)
                sample = str.split(file,'_')[0]
                res_dir = Folder + con + '_ML_' + str(sample) + '.csv'
                toapp_dir = Folder + con + '_ML_' + file + '.csv'
                if os.path.exists(res_dir):
                    Tomod = pd.read_csv(res_dir) 
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    #Count=Count.drop_duplicates()
                    #print(Count)
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Tomod = Tomod.append(Toappend)
                    Tomod.to_csv(res_dir,index = False,header=True)
                    os.remove(toapp_dir)
                else:
                    Toappend = pd.read_csv(toapp_dir)
                    Toappend['bin'] = [((x-1)//400)*400+200  for x in Toappend['pos']]
                    Count = Toappend.groupby(['chrom','bin','strand']).size().reset_index(name='counts')
                    Toappend=Toappend.merge(Count, on=['chrom','bin','strand'])
                    #print(Toappend)
                    conditions = [
                        (Toappend['counts'] > 4),
                        (Toappend['counts'] < 5) 
                    ]
                    # create a list of the values we want to assign for each condition
                    values = [Toappend['ML'], np.nan]

                    # create a new column and use np.select to assign values to it using our lists as arguments
                    Toappend['ML'] = np.select(conditions, values)
                    Toappend=Toappend.drop(columns=['counts','pos'])
                    #Toappend=Toappend.dropna(axis = 0, thresh=4, inplace = True)
                    Toappend=Toappend.groupby(['chrom','bin','strand']).agg({'ML': 'mean'}).reset_index()
                    Toappend.to_csv(res_dir,index = False,header=True)
                    os.remove(toapp_dir)
                    
        logm("Merging MeH between samples for CHH.")            
        # merge MeH between samples
        for sample in bam_list: 
            tomerge_dir = Folder +  con + '_' + str(sample) + '.csv' 
            res_dir = Folder +  con + '_' + 'Results.csv'
            if os.path.exists(res_dir):
                Result = pd.read_csv(res_dir)
                Tomerge = pd.read_csv(tomerge_dir)
                #Tomerge = Tomerge.drop(columns=['dis','ML','depth'])
                Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                Tomerge = Tomerge.rename(columns={'MeH': sample})
                Result = Result.merge(Tomerge, on=['chrom','bin','strand'])
                Result.dropna(axis = 0, thresh=4, inplace = True) 
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
            else:
                Result = pd.read_csv(tomerge_dir)
                Result.head()
                #Result = Result.drop(columns=['dis','ML','depth'])
                Result.dropna(axis = 0, thresh=4, inplace = True)
                Result = Result.rename(columns={'MeH': sample})
                Result.to_csv(Folder + con + '_' +'Results.csv',index = False,header=True)
                os.remove(tomerge_dir)
        logm("Merging ML between samples for CHH.")    
        # merge ML between samples
        if args.mlv:
            for sample in bam_list: 
                tomerge_dir = Folder +  con + '_ML_' + str(sample) + '.csv' 
                res_dir = Folder +  con + '_ML_' + 'Results.csv'
                if os.path.exists(res_dir):
                    Result = pd.read_csv(res_dir) 
                    Tomerge = pd.read_csv(tomerge_dir)
                    Tomerge.dropna(axis = 0, thresh=4, inplace = True)
                    Tomerge = Tomerge.rename(columns={'ML': sample})
                    Result=Result.merge(Tomerge, on=['chrom','bin','strand'])
                    Result.dropna(axis = 0, thresh=4, inplace = True) 
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
                else:
                    Result = pd.read_csv(tomerge_dir)
                    Result = Result.rename(columns={'ML': sample})
                    #Result = Result.drop(columns=['counts','pos','depth','dis'])
                    Result.dropna(axis = 0, thresh=4, inplace = True)
                    Result.to_csv(res_dir,index = False,header=True)
                    os.remove(tomerge_dir)
        #Result.to_csv(res_dir ,index = False,header=True)
        print("All done.",len(bam_list),"bam files processed and merged for CHH.")    
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CHH.")
        
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)

    topp=topp.groupby(['context','sample']).agg({'context_coverage': 'sum', 'coverage': 'sum'}).reset_index()
    #print('after groupby',topp) 
    
    for filename in spbam_list:
        file = Folder + filename + '.bam'
        fileind = Folder + filename + '.bam.bai'
        os.remove(file)
        os.remove(fileind)
        
    
    end = time.time()
    
    for i in range(topp.shape[0]):
        #print('i = ',i)
        print('Sample', topp.iloc[i,1],'has coverage',topp.iloc[i,2],'for context',topp.iloc[i,0],'out of data coverage',topp.iloc[i,3])
        logm('Sample '+str(topp.iloc[i,1])+' has coverage '+str(topp.iloc[i,2])+' for context '+str(topp.iloc[i,0])+' out of data coverage '+str(topp.iloc[i,3])+ '.')


# FINAL FINAL
# /MH/test
# python3 testfull.py -w 4 -c 4 --CG --CHG --opt --mlv
# python3 finalfinal.py -w 4 -c 80 --CG --opt --mlv
# MH/testsp/hg19
# python testfull.py -w 4 -c 8 --CG







