import os
import re
import numpy
import math
import random
import matplotlib.pyplot as plt

####P+先验概率构建一个16*8的矩阵A,C,G,T.  AA=0 AC=1 AG=2 AT=3  CA=4 ... TT=15
####P-后验概率构建一个16*8的矩阵A,C,G,T.  AA=0 AC=1 AG=2 AT=3  CA=4 ... TT=15

count_test=0
def fuzhi(S):
    if S=="A" :
        a=0
        return a
    if S=="C":
        a=1
        return a
    if S=="G":
        a=2
        return a
    if S=="T":
        a=3
        return a
    return a

def fenbu(donor):
    count=numpy.zeros((9,4))
    Pp=numpy.zeros((8,16))
    Pfu=numpy.zeros((8,16))
    for i in range(len(donor)):
        for m in range(len(donor[i])):
            if donor[i][m]=="A":
                count[m][0]+=1
            if donor[i][m]=="C":
                count[m][1]+=1
            if donor[i][m]=="G":
                count[m][2]+=1
            if donor[i][m]=="T":
                count[m][3]+=1
            if m<8:
                b=4*fuzhi(donor[i][m])+fuzhi(donor[i][m+1])
                Pp[m][b] +=1
    for i in range(len(Pp)):
        for j in range(len(Pp[i])):
            if count[i][j//4]/max(count[i][0],count[i][1],count[i][2],count[i][3])<0.1 and i !=2:
                Pfu[i][j]=0
               
            else:
                Pfu[i][j]=(Pp[i][j]/count[i][j//4]).round(4)   
    return count,Pp,Pfu

def jisuanS(P1,p0,P2,p,s):
    Sx=0
    for n in range(len(s)):
        if s[n]!="A" and s[n]!="C"and s[n]!="G" and s[n]!="T" :
            Sx=1000
            return Sx
    for i in range(len(s)):
        b=fuzhi(s[i-1])*4+fuzhi(s[i])
        if i==0:
            Sx=math.log(p0[fuzhi(s[0])]/p[fuzhi(s[0])])+Sx
        elif  P2[i-1][b]==0 or P1[i-1][b]==0:
            Sx=Sx+math.log(1)
        else:                    
            Sx=math.log(P1[i-1][b]/P2[i-1][b])+Sx
    return Sx   

def panduan(s):
    a=0
    for i in range(len(s)):
        if s[i]!="A" and s[i]!="C" and s[i]!="G" and s[i]!="T" :
            a=1
    return a   

j=0
jj=0
Donor={}
d_train={}
TP=0
FN=0
FP=0
count_test=0
root_dir=r"S:\Data\Training Set"
if True:
    for file in os.listdir(root_dir):
        file_name=root_dir+ "\\" + file
        filein=open(file_name,"r")
        wei=[]
        Ss=[]
        xulie=str()
        
        for line in filein:
            if line[0:3]=="CDS":
                weidian=re.sub("\D", " ", line)
                #print(line)
                weidian=weidian.split()
                
                for i in range (len(weidian)):
                    if i%2 != 0:
                        wei.append(weidian[i])
            if line[0]=="a" or line[0]== "t" or line[0]== "c" or line[0]== "g":
                xulie=xulie+line[0:-1]
        
        for i in range (len(wei)):
            for m in range(5):
                aaaa=random.randint(1,len(xulie)-9)
                Ss=xulie[aaaa:aaaa+9].upper()
                aa=panduan(Ss)
                if aa==0:
                    d_train[jj]=Ss
                    jj=jj+1
            Ss=xulie[eval(wei[i])-3:eval(wei[i])+6].upper()
            if  Ss[1]!="N" and len(Ss)==9 :
                Donor[j]=Ss
                j=j+1
    
                
            ####读出donor位点的序列

count,count_zong,Pzheng=fenbu(Donor)
count_fu,count_zong,Pfu=fenbu(d_train)
Pzheng_1=count[0]/j
Pfu_1=count_fu[0]/jj

Sx=0
for i in range(len(Donor)):
    Sx=jisuanS(Pzheng,Pzheng_1,Pfu,Pfu_1,Donor[i])+Sx
#print(Sx/j)
print(Sx/j)
aaaaa=Sx/j
print(Pzheng,"\n",Pzheng_1,"\n",Pfu,"\n",Pfu_1)
print("##############")
jjj=0
jjjj=0
Donor_test={}
d_test={}
root_dir=r"S:\Data\Testing Set"
if True:
    for file in os.listdir(root_dir):
        file_name=root_dir+ "\\" + file
        filein=open(file_name,"r")
        xulie_test=str()
        wei=[]
        for line in filein:
            if line[0]==">" and len(line)> 25:
                weidian1=str(re.findall(r'[(](.*?)[)]', line))
                weidian=re.sub("\D", " ", weidian1)
                weidian=weidian.split()
                for i in range (len(weidian)):
                    if i%2 != 0:
                        #print(weidian[i])
                        wei.append(weidian[i])
                #####找出donor位点
            if line[0]=="A" or line[0]== "T" or line[0]== "C" or line[0]== "G":
                xulie_test=xulie_test+line[0:-1]
        for i in range (len(wei)):
            for m in range(3):
                aaaa=random.randint(1,len(xulie_test)-9)
                Ss=xulie_test[aaaa:aaaa+9].upper()
                aa=panduan(Ss)
                if aa==0:
                    d_test[jjj]=Ss
                    jjj=jjj+1
            Ss=xulie_test[eval(wei[i])-3:eval(wei[i])+6].upper()
            aa=panduan(Ss)
            if  aa==0 and eval(wei[i])<(len(xulie_test)-6):
                    Donor_test[jjjj]=Ss
                    jjjj=jjjj+1
                  
####Sn=TP/(TP+FN) TP表示预测的真实的Donor剪接位点,FN表示未能预测到的但仍是真实的Donor剪接位点，
####Sp=TP/(TP+FP) FP表示预测的假的Donor位点
C=0
print(jjjj,jjj)
XSn=numpy.zeros(100)
XSp=numpy.zeros(100)
aa=100
for C in range(11,111):
    TP=0
    FP=0
    FN=0
    for i in range (jjj):
        Cfen=jisuanS(Pzheng,Pzheng_1,Pfu,Pfu_1,d_test[i])
        if abs(Cfen-aaaaa)<C/18:
       # if math.log(abs(Cfen-aaaaa)/aaaaa)<C/50:
            FP=FP+1
        
    for i in range(jjjj):
        #print(Donor_test[i])
        Cfen=jisuanS(Pzheng,Pzheng_1,Pfu,Pfu_1,Donor_test[i])
        if abs(Cfen-aaaaa)<C/18:
       # if  math.log(abs(Cfen-aaaaa)/aaaaa)<C/50:
            TP=TP+1
    FN=jjjj-TP
    Sn=TP/(TP+FN)
    Sp=TP/(TP+FP)
    aa=aa-1
    XSn[aa]=Sn
    XSp[aa]=Sp             

print(XSn,XSp)
# 绘制余弦曲线，使用蓝色的、连续的、宽度为 1 （像素）的线条
plt.plot(XSn[0:-1], XSp[0:-1],"-")
plt.xlabel("Sn")
plt.ylabel("Sp")
plt.show()