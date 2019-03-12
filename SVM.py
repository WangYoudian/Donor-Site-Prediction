# coding:utf-8
import numpy
import os
import re
import random
from sklearn import svm
import matplotlib.pyplot as plt
def trans(s):
    aaa=[]
    it={"A": [0,0,0,1] ,"G":[0,0,1,0],"C":[0,1,0,0],"T":[1,0,0,0]}
    for i in range(len(s)):
        aaa=aaa+it[s[i]]   
    return aaa
def panduan(s):
    a=0
    for i in range(len(s)):
        if s[i]!="A" and s[i]!="C" and s[i]!="G" and s[i]!="T" :
            a=1
    return a

Donor_train={}
d_train={}
Donor_test={}
d_test={}
TP=0
FN=0
FP=0
j=0
jj=0
jjj=0
jjjj=0
root_dir=r"S:\Data\Training Set"
#with open("Training set.txt","w") as train:
if True:
    for file in os.listdir(root_dir):
        file_name=root_dir+ "\\" + file
        filein=open(file_name,"r")
        wei=[]
        xulie=str()
        ####批量读入文件夹下的txt文件
        for line in filein:
            if line[0:3]=="CDS":
            	##
                weidian=re.sub("\D", " ", line)
                weidian=weidian.split()
                for i in range (len(weidian)):
                    if i%2 != 0:
                        wei.append(weidian[i])
            if line[0]=="a" or line[0]== "t" or line[0]== "c" or line[0]== "g":
                xulie=xulie+line[0:-1]
        print(wei)
        for i in range (len(wei)):
            for m in range(3):
                aaaa=random.randint(1,len(xulie)-9)
                Ss=xulie[aaaa:aaaa+9].upper()
                aa=panduan(Ss)
                if aa==0:
                    d_train[jj]=Ss
                    jj=jj+1
            Ss=xulie[eval(wei[i])-3:eval(wei[i])+6].upper()
            if  Ss[1]!="N" and len(Ss)==9 :
                    Donor_train[j]=Ss
                    j=j+1
			#print(Donor_train[j-1],j,len(Donor_train))
			#aaa=trans(d_train[1])
			#print(aaa,jj)            

# Parsing the data       
Dtrain=numpy.zeros(shape=(j+jj,36))
#dtrain=numpy.zeros(shape=(jj,36))
y=numpy.zeros(j+jj)
for i in range(len(Donor_train)):
    Dtrain[i]=trans(Donor_train[i])
    y[i]=1
for i in range(len(d_train)):
    Dtrain[i+j]=trans(d_train[i])
    
print(Dtrain,d_train)
#clf = SVC(C=0.8, kernel='rbf', gamma=20, decision_function_shape='ovo')
#clf.fit(Dtrain[1],dtrain[1])

#X = [[0, 0], [1, 1],[1,0],[1,1],[2,2],[1,2]]
#yy = [0, 1,1,1,1,1]
#clf = svm.SVR(kernel='poly',C=1,degree=2,gamma=0.5)
clf = svm.SVR(kernel='rbf',gamma=0.5)
clf.fit(Dtrain, y)
result=clf.predict(Dtrain[j-50:j+50])
print(j,jj)
print(result)

####读出donor位点的序列
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
                        wei.append(weidian[i])
                
            if line[0]=="A" or line[0]== "T" or line[0]== "C" or line[0]== "G":
                xulie_test=xulie_test+line[0:-1]
        for i in range (len(wei)):
            for m in range(6):
                aaaa=random.randint(1,len(xulie_test)-9)
                Ss=xulie_test[aaaa:aaaa+9].upper()
                aa=panduan(Ss)
                if aa==0:
                    d_test[jjj]=Ss
                    jjj=jjj+1
            Ss=xulie_test[eval(wei[i])-3:eval(wei[i])+6].upper()
            aa=panduan(Ss)
            if  aa==0 and eval(wei[i])<(len(xulie_test)-6) :
                    Donor_test[jjjj]=Ss
                    jjjj=jjjj+1
####Sn=TP/(TP+FN) TP表示预测的真实的Donor剪接位点,FN表示未能预测到的但仍是真实的Donor剪接位点，
####Sp=TP/(TP+FP) FP表示预测的假的Donor位点

XSp=numpy.zeros(50)
XSn=numpy.zeros(50)
aa=0
print(j,jjjj)
for C in range(50):
    TP=0
    FP=0
    FN=0
    for i in range (jjj):
        Cfen=clf.predict([trans(d_test[i])])
        if Cfen>C/51:
            FP=FP+1
        
    for i in range( jjjj):
        Cfen=clf.predict([trans(Donor_test[i])])
        if Cfen>C/51:
            TP=TP+1
    FN=jjjj-TP
    Sn=TP/(TP+FN)
    Sp=TP/(TP+FP)
    XSp[aa]=Sp
    XSn[aa]=Sn
    aa+=1

print(XSn,XSp)
plt.plot(XSn, XSp, color="blue", linewidth=1.0, linestyle="-")
plt.xlabel("Sn")
plt.ylabel("Sp")
plt.show()
