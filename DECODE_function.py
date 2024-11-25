#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os, re, csv, time,struct,time,subprocess,base64,svgwrite,warnings, copy
from lttb import largest_triangle_three_buckets
import numpy as np
import pandas as pd
from scipy import  misc
from scipy import ndimage as ndi
from scipy import signal
from scipy import stats
import scipy.stats as st
import matplotlib.pyplot as plt
from PIL import Image
warnings.filterwarnings('ignore')

LOGO_start=0
LOGO_length =12

AA = "ARNDCQEGHILKMFPSTWYV*"
cname=["ProteinID","pos","AA","read","distance"]

FPaa=r'C:\Users\DECODE\Desktop\AA'
AAlist_logo="ARNDCQEGHILKMFPSTWYV"

imlist=[]
for i in range(len(AAlist_logo)):
    imlist.append(Image.open(FPaa+"\\"+AAlist_logo[i]+".png"))
brank =Image.open(FPaa+"\\Brank.png")

AAlist=[a for a in AA]
cname = cname+AAlist


dt_AWAS = np.dtype([
    ('ProteinID', 'i4'), ('pos', 'i4'), ('AA', 'i4'), ('read', 'f4'), ('distance', 'f4'),
    ('A', 'f4'),('R', 'f4'), ('N', 'f4'),('D', 'f4'), ('C', 'f4'),('Q', 'f4'), ('E', 'f4'),('G', 'f4'), ('H', 'f4'),('I', 'f4'), 
    ('L', 'f4'),('K', 'f4'), ('M', 'f4'),('F', 'f4'), ('P', 'f4'), ('S', 'f4'),('T', 'f4'), ('W', 'f4'),('Y', 'f4'), ('V', 'f4'),('*', 'f4')
])

dt_AWAS_Protein_rank = np.dtype([
     ('barcode', 'i4'), ('ID', 'i4'), ('pos', 'i4'), ('rank', 'i4'), ('read', 'f4'), ('distance', 'f4'),])



def thinout_array(array ,width):
    array2 = ndi.maximum_filter(array, size=width, mode='nearest')
    array2 = array2[::width]
    return array2

def thinout_index(array,width):
    array = array[::width]
    return array

def SortAWASPlot100(df,rp,wp,fname,targetID,dswidth, plotdata):
    
    x=[i for i in range(100)]
    y=[df.iloc[i][plotdata] for i in range(100)]
    cmap =[df.iloc[i]["color"] for i in range(100)]
    return x,y,cmap

def SortAWASPlot1000(df,rp,wp,fname,targetID,dswidth, plotdata):
    
    x=[i for i in range(1000)]
    y=[df.iloc[i][plotdata] for i in range(1000)]
    cmap =[df.iloc[i]["color"] for i in range(1000)]
    return x,y,cmap

def ReadTotalRead(fp):
    totalread=0
    with open(fp, 'r', encoding='UTF-8') as f:
        for line in f:
            line=line.split(",")
            totalread+=int(line[1])
    print("total read:"+str(totalread)+" reads")
    return totalread

def MakeAWASPlot(df,rp,wp,fname,targetID,dswidth, readlimit=0.001):
    df0=df.loc[0:df.index[df['ProteinID']==targetID].tolist()[0]]
    dft=df.loc[df.index[df['ProteinID']==targetID].tolist()[0]:df.index[df['ProteinID']==targetID].tolist()[len(df.index[df['ProteinID']==targetID].tolist())-1]]
    df1=df.loc[df.index[df['ProteinID']==targetID].tolist()[len(df.index[df['ProteinID']==targetID].tolist())-1]:len(df)-1]
  
    nmp=list(np.array(df0["data"]))
    tdata0 = [(i,nmp[i]) for i in  range(len(df0))]
    sampled_dataset0 = largest_triangle_three_buckets(tdata0, int(len(tdata0)/dswidth))
    cmap0 = ["gray"] * len(sampled_dataset0)

    nmp=list(np.array(dft["data"]))
    tdatat = [(i,nmp[i]) for i in  range(len(dft))]
    if(dswidth>len(tdatat)):
        print ("Too shot target length:"+(len(tdatat)))
        return -1
    sampled_datasett = largest_triangle_three_buckets(tdatat, int(len(tdatat)/dswidth))
    cmapt = ["red"] * len(sampled_datasett)

    nmp=list(np.array(df1["data"]))
    tdata1 = [(i,nmp[i]) for i in  range(len(df1))]
    sampled_dataset1 = largest_triangle_three_buckets(tdata1, int(len(tdata1)/dswidth))
    cmap1 = ["gray"] * len(sampled_dataset1 )
   
    sampled_dataset = sampled_dataset0+sampled_datasett+sampled_dataset1
    cmap = cmap0+cmapt+cmap1
    x=[i for i in range(len(sampled_dataset))]
    y=[i[1] for i in sampled_dataset]

    return x,y,cmap
    
    
def MakeAWASPlot_woTarget(df,rp,wp,fname,targetID,dswidth):
    nmp=list(np.array(df["data"]))
    tdata = [(i,nmp[i]) for i in  range(len(df))]
    sampled_dataset = largest_triangle_three_buckets(tdata, int(len(tdata)/dswidth))
    cmap = ["gray"] * len(sampled_dataset)
    x=[i for i in range(len(sampled_dataset))]
    y=[i[1] for i in sampled_dataset]

    return x,y,cmap
    
def MakePlot_Target(rp,wp,fname,targetID,dswidth, plotdata, sort, mean,figinfo=[14,10,8,10,0,0,0,0], rlimit=1):
    #plotdata="read"."distance","readxdistance"
    #figinfo=[fontsize1,width,heigh,dot-size,xlim,ylim,num-of-xticks,num-of-yticks ]
    if(len(figinfo)!=8):
        print("Error arg of figinfo list ")
        print("figinfo=[fontsize1,width,heigh,dot-size,xlim,ylim,num-of-xticks,num-of-yticks]")
        return -1
    fontsize1 = figinfo[0]
    wpd=wp+"\\"+plotdata
    print("targetID:" + str(targetID))
    d =  np.fromfile(rp+"\\"+fname, dtype=dt_AWAS)
    df=pd.DataFrame(d)
    df.columns=cname
    df=df.drop(0)
    df.fillna("0A")
    df['ProteinID']=df['ProteinID'].astype(int)
    df['color'] = "gray"
    if(targetID !=-1):
        df["color"]=df['color'].mask(df['ProteinID'] == targetID,"red")
    df.fillna(0)
    df["pos"]=df["pos"].astype(float)
    df["pos"]=df["pos"].astype(int)
    df["read"]=df["read"].astype(float)
    df["distance"]=df["distance"].astype(float)
   
    df.loc[df["read"]<rlimit, ("read","distance")]=0#特定のread数以下を0に   
    v = (np.ones(5))/5.0  
    
    if(plotdata =="readxdistance"):
        nmp1=np.array(df["read"])
        nmp2=np.array(df["distance"])
        nmp= nmp1* nmp2
        df["data"]=  np.convolve(nmp, v, mode='same')
    else:
        if(mean==1):
            nmp=np.array(df[plotdata])
            df["data"]=  np.convolve(nmp, v, mode='same')
            wpd=wpd+"_mean"
        else:
            df["data"]= df[plotdata]

    if(sort==1):
        df=df.sort_values("data", ascending=False)
        wpd=wpd+"_sort"
     
    if not os.path.exists(wpd+"\\plot"):
        os.makedirs(wpd+"\\plot")
    if not os.path.exists(wpd+"\\csv"):
        os.makedirs(wpd+"\\csv")

    if(sort==1):
        result = SortAWASPlot100(df,rp,wp,fname,targetID,dswidth)
    else:
        if(targetID !=-1):
            result = MakeAWASPlot(df,rp,wp,fname,targetID,dswidth)
        else:
            result = MakeAWASPlot_woTarget(df,rp,wp,fname,targetID,dswidth)
    if(result ==-1):
        return -1
    plt.figure(figsize=(figinfo[1],figinfo[2]))
    plt.title(fname+"AWAS plot",fontsize= fontsize1)
    plt.xlabel("protein(A.A.)")
    plt.ylabel(plotdata)
    plt.scatter(result[0], result[1], s=figinfo[3],c=result[2])   
    if(figinfo[4] !=0):
        plt.xlim(0, figinfo[4] )
    if(figinfo[6] !=0):
        if(figinfo[4] !=0):
            plt.xticks(np.arange(0, figinfo[4],int(figinfo[4]/figinfo[6])))
        else:
            plt.xticks(np.arange(0, len(result[0]),int(len(result[0])/figinfo[6])))
        
    if(figinfo[5] !=0):
        plt.ylim(0, figinfo[5] )
        
    if(figinfo[7] !=0):
        if(figinfo[5] !=0):
            plt.yticks(np.arange(0, figinfo[5],int(figinfo[5]/figinfo[7])))
        else:
            maxval=max(result[0])
           # print("maxval: "+str(maxval))
            plt.yticks(np.arange(0, maxval,int(maxval/figinfo[7])))
        
    plt.savefig(wpd+"\\plot\\"+fname+"_AWAS.svg")
    plt.savefig(wpd+"\\plot\\"+fname+"_AWAS.png")
    plt.show()
    
    
    plt.figure(figsize=(figinfo[1],figinfo[2]))
    plt.title(fname+"AWAS plot",fontsize= fontsize1)
    plt.xlabel("protein(A.A.)")
    plt.ylabel(plotdata)
    plt.scatter(result[0], result[1], s=figinfo[3],c=result[2])   
    if(figinfo[4] !=0):
        plt.xlim(0, figinfo[4] )
    if(figinfo[6] !=0):
        if(figinfo[4] !=0):
            plt.xticks(np.arange(0, figinfo[4],int(figinfo[4]/figinfo[6])))
        else:
            plt.xticks(np.arange(0, len(result[0]),int(len(result[0])/figinfo[6])))
        
    if(figinfo[5] !=0):
        plt.ylim(0, figinfo[5] )
        
    if(figinfo[7] !=0):
        if(figinfo[5] !=0):
            plt.yticks(np.arange(0, figinfo[5],int(figinfo[5]/figinfo[7])))
        else:
            maxval=max(result[0])
           # print("maxval: "+str(maxval))
            plt.yticks(np.arange(0, maxval,int(maxval/figinfo[7])))
            
    plt.yscale("log")
    plt.savefig(wpd+"\\plot\\"+fname+"_AWAS_log.svg")
    plt.savefig(wpd+"\\plot\\"+fname+"_AWAS_log.png")
    plt.show()
   
    dft=df[df['ProteinID']==targetID] 
    dft=dft.reset_index()
    
    #exdf=dft.filter(items=['ProteinID', "pos","data","read","distance"], axis=1)
    #exdf.to_csv(wpd+"\\csv\\"+fname+".csv");
    dft.to_csv(wpd+"\\csv\\"+fname+".csv");
    
    maxid = signal.argrelmax(np.array(dft["data"]), order=12) 
    peakdf=dft.loc[maxid].sort_values("data", ascending=False)
    #print(peakdf)
    lim=0

    for i in range(len(peakdf)):
        if(peakdf.iloc[i]["read"]<rlimit):
            continue
        lim+=1
        
        if( lim>10):
            break
        print(peakdf.iloc[i]["index"])
        maxpos=peakdf.iloc[i]["index"]
        LOGO_start=peakdf.iloc[i]["index"]
        if(df.shape[0]<=maxpos+LOGO_length):
            LOGO_start=df.shape[0]-LOGO_length-1
        AAratio=np.zeros((LOGO_length,20))
        for a in range(20):
            for p in range(LOGO_length):

                AAratio[p][a]=df[AAlist_logo[a]].iloc[LOGO_start+p]
        AAratio=AAratio*256

        LOGO=imlist[9].resize((25,256)) # make left scale using I
        #print(AAratio.shape)
        for p in range(AAratio.shape[0]):
            ResultAA=brank.resize((256,1))
            for n in range(AAratio.shape[1]):
                if(AAratio[p][n]>256*20/100):#20% 以上を表示
                #if(AAratio[p][n]>25):
                    ResultAA=get_concat_v(ResultAA, imlist[n].resize((256,int(AAratio[p][n]))))
            LOGO=get_concat_h(LOGO,ResultAA)
        #print(num+"_target:"+str(targetID)+"_"+str(LOGO_start))
        plt.imshow(np.array(LOGO))
        if not os.path.exists(wpd+"\\Logo\\"+fname[0:-4]):
            os.makedirs(wpd+"\\Logo\\"+fname[0:-4])
        LOGO.save(wpd+"\\Logo\\"+fname[0:-4]+"\\target_"+str(targetID)+"_Rank"+str(i)+"_"+str(peakdf.iloc[i]["pos"])+"_Score"+str(peakdf.iloc[i]["data"])+".png", quality=100)
    return dft,peakdf
        
def MakeLogo(rp,wp,fname,targetID, plotdata,logomin=20):
    fontsize1 = 14
    print("targetID:" + str(targetID))
    d =  np.fromfile(rp+"\\"+fname, dtype=dt_AWAS)
    df=pd.DataFrame(d)
    df.columns=cname
    df=df.drop(0)
    df.fillna("0A")
    df['ProteinID']=df['ProteinID'].astype(int)

    dft=df[df['ProteinID']==targetID]
    dft.fillna(0)
    dft["pos"]=dft["pos"].astype(float)
    dft["pos"]=dft["pos"].astype(int)
    dft["read"]=dft["read"].astype(float)
    dft["distance"]=dft["distance"].astype(float)
    
    
    wpd=wp+"\\Logo"

    dft=dft.reset_index()
    
    maxpos=dft[plotdata].idxmax(axis = 0)
    print("maxpos: "+str(maxpos))
  
    LOGO_start=dft.iloc[maxpos]["index"]
    if(df.shape[0]<=maxpos+LOGO_length):
        LOGO_start=df.shape[0]-LOGO_length-1
    AAratio=np.zeros((LOGO_length,20))
    for a in range(20):
        for p in range(LOGO_length):
            #print(str(a)+"_"+str(p)+"_"+str(LOGO_start+p))
            AAratio[p][a]=df[AAlist_logo[a]].iloc[int(LOGO_start+p)]
    AAratio=AAratio*256

    LOGO=imlist[9].resize((25,256)) # make left scale using I
    #  print(AAratio.shape)
    for p in range(AAratio.shape[0]):
        ResultAA=brank.resize((256,1))
        for n in range(AAratio.shape[1]):
            if(AAratio[p][n]>256*logomin/100):#logomin% 以上を表示
            #if(AAratio[p][n]>25):
                ResultAA=get_concat_v(ResultAA, imlist[n].resize((256,int(AAratio[p][n]))))
        LOGO=get_concat_h(LOGO,ResultAA)
    print("target:"+str(targetID)+"_"+str(LOGO_start))
    plt.imshow(np.array(LOGO))
    
    if not os.path.exists(wpd):
        os.makedirs(wpd)
  
     
    LOGO.save(wpd+"\\"+fname[0:-4]+"_target_"+str(targetID)+"_"+str(maxpos)+".png", quality=100)
  
    return dft
def matching(s,t):
    if(len(s)<len(t)):
        return -1
    count=[0 for i in  range (len(s))]
    for p in range(len(s)-len(t)):
        for a in range(len(t)):
            if(s[p+a]==t[a]):
                count[p]+=1
   # print(count)
    if(max(count) >3):    
        return count.index(max(count))
    else:
        return -1
   
    return -1            
def MakeLogo_Tag(Tag,FPr,fname,mode=0,logomin=20,wpd="./temp"):
    LOGO_length=len(Tag)
    total_read=0

    csv_file = open(FPr+"\\"+fname, "r", encoding="ms932", errors="", newline="" )
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)


    AAratio=np.zeros((LOGO_length,len(AAlist_logo)))
    for row in f:
        #print(row)
        total_read +=int(row[1])
        pstart=matching(row[0],Tag)
        """matchObj = re.search(Tag, row[0])
        if(len(row[0])<=LOGO_length):
            continue
        if (matchObj):
            #print (matchObj.group()) # マッチした文字列： abc
           # print (matchObj.start()) # マッチした文字列の開始位置： 3
           # print (matchObj.span()) # マッチした文字列の開始位置と終了位置： (3, 6)
            pstart=int(matchObj.start())
          """
        #print(pstart)
        if(pstart==-1):
            continue
           
        if(mode==0):
            for p in range(LOGO_length-pstart):
                for a in range(len(AAlist_logo)):
                    #print(row[0][p+pstart]+":"+a)
                    if(row[0][p+pstart]==AAlist_logo[a]):
                        AAratio[p][a]+=int(row[1])
                        #print(row)
                        break
        elif(mode==1):
            for p in range(LOGO_length-pstart):
                for a in range(len(AAlist_logo)):
                    #print(row[0][p+pstart]+":"+a)
                    if(row[0][p+pstart]==AAlist_logo[a]):
                        AAratio[p][a]+=1
                        #print(row)
                        break                           

        #break

    for a in range(len(AAratio)):
        norm=0
        for n in range(len(AAratio[a])):
            #print(AAratio[a][n])
            norm+=AAratio[a][n]
        #print(norm)
        AAratio[a]=AAratio[a]/norm
    AAratio=AAratio*256

    LOGO=imlist[9].resize((25,256)) # make left scale using I
    #  print(AAratio.shape)
    for p in range(AAratio.shape[0]):
        ResultAA=brank.resize((256,1))
        for n in range(AAratio.shape[1]):
            if(AAratio[p][n]>256*logomin/100):#min% 以上を表示
            #if(AAratio[p][n]>25):
                ResultAA=get_concat_v(ResultAA, imlist[n].resize((256,int(AAratio[p][n]))))
        LOGO=get_concat_h(LOGO,ResultAA)
    plt.imshow(np.array(LOGO))        

    if not os.path.exists(wpd):
        os.makedirs(wpd)

    df = pd.DataFrame(AAratio) 
    df.to_csv(wpd+"\\"+fname[0:-4]+"_Tag_"+Tag.replace('.', 'x')+".csv")
    LOGO.save(wpd+"\\"+fname[0:-4]+"_Tag_"+Tag.replace('.', 'x')+".png", quality=100)
    return AAratio
    
def Extract_target_protein(rp,targetID):
    print("import: "+rp)
    d =  np.fromfile(rp, dtype=dt_AWAS)
    df=pd.DataFrame(d)
    df.columns=cname

    df=df.drop(0)
    df.fillna("0A")
    df['ProteinID']=df['ProteinID'].astype(int)
    df=df.fillna(0)
    df_target=df[df['ProteinID']==targetID]

    return df_target

def Extract_target_protein_save(rp,wp,fname,targetID):
    d =  np.fromfile(rp+"\\"+fname, dtype=dt_AWAS)
    df=pd.DataFrame(d)
    df.columns=cname

    df=df.drop(0)
    df.fillna("0A")
    df['ProteinID']=df['ProteinID'].astype(int)
    df=df.fillna(0)
    df_target=df[df['ProteinID']==targetID]
    
    if not os.path.exists(wp+"\\"+fname):
        os.makedirs(wp+"\\"+fname) 
    df_target.to_csv(wp+"\\"+fname+"\\"+str(targetID)+".csv")
    print("export:"+wp+"\\"+fname+"\\"+str(targetID)+".csv")
    return df_target
    
    
def minmax_norm(df):
    df2=df.copy()
    for ff in df:
        #print(ff+str(df[ff].min())+str( df[ff].max()))
        df2[ff]=(df[ff] - df[ff].min()) / ( df[ff].max() - df[ff].min())
    return df2

def get_concat_h(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, im1.height),"white")
    dst.paste(im1, (0, dst.height-im1.height))
    dst.paste(im2, (im1.width, dst.height-im2.height))
    return dst

def get_concat_v(im1, im2):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height),"white")
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst 

def Fastq2Pep(Fastq2Pep_exe,FP_fastq,FP_Pep,Codon,Barcode):
    cmd = Fastq2Pep_exe+" "+FP_fastq+" "+FP_Pep+" "+Codon+" "+Barcode
    awas=subprocess.Popen(cmd.split(), shell=True)
    time_sta = time.time()
    while(awas.poll()==None):
        time.sleep(60) 
        time_end = time.time()
        tim = time_end- time_sta
        tim=tim/60
        if(int(tim%10)==0):
            print(str(tim)+" min")

def AWASanalysis(AWAS_exe_cuda,FP_Pep,AWAS_OUT,P_MAP):
    peplist_d = [d for d in os.listdir(FP_Pep)]
    for d in peplist_d:
        Pep_IN=FP_Pep+"\\"+d+"\\Peptide"
        AWAS_OUT_t=AWAS_OUT+"\\"+d
        cmd = AWAS_exe_cuda\
        +r" --pep "+Pep_IN \
        +r" --export "+AWAS_OUT_t \
        +r" --P_map "+P_MAP\
        +" --Calc read"\
        +" --readlimit 0"\
        +" --Distance_function 3"\
        +" --TH 0.001"\
        +" --Start 0 "
        cmd=cmd.split()
        time_sta = time.time()
        while(1):
            awas=subprocess.Popen(cmd, shell=True)
            while(awas.poll()==None):
                time.sleep(60) 
                time_end = time.time()
                tim = time_end- time_sta
                tim=tim/60
                if(int(tim%10)==0):
                    print(str(tim)+" min")
            data=[]
            with open(r'AWAS_Running.log', 'r', encoding='UTF-8') as f:
                data = f.read()
                data=data.split("\n")
                print(data)
            if(int(data[0])==(int(data[1])-1)):
                break
            else:
                cmd[16]=data[0]

# In[ ]:




