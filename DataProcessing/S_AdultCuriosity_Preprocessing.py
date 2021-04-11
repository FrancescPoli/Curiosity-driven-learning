# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:32:15 2020

@author: Francesco
"""
import os
import numpy as np
import pandas as pd
import glob
import math
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

# Set directory
os.chdir(r"U:\Desktop\S_AdultCuriosity\DataProcessing\data\online_game")
#Import all csv files in the directory as one dataframe
l = [pd.read_csv(filename) for filename in glob.glob("*.csv")]
df = pd.concat(l, axis=0)
del(l)

# Remove rows that are not needed
df.dropna(subset = ["mouse.time"], inplace=True)
# Remove cols that are not needed
df=df.drop(['key_resp.keys', 'key_resp.rt', 'mouse_2.x', 'mouse_2.y', 'mouse_2.rightButton',
         'mouse_2.leftButton', 'mouse_2.midButton', 'mouse.leftButton',
         'mouse.midButton', 'mouse.rightButton', 'mouse.clicked_pos', 'sequenceloop.thisRepN',
         'sequenceloop.thisTrialN', 'sequenceloop.thisN',
         'sequenceloop.thisIndex', 'sequenceloop.ran', 'sequenceloop.order',
         'change_character.thisRepN', 'change_character.thisTrialN',
         'change_character.thisN', 'change_character.thisIndex',
         'change_character.ran', 'change_character.order', 'date', 'expName', 'psychopyVersion', 'OS', 'frameRate',
         'target_loop.thisRepN', 'target_loop.thisTrialN', 'target_loop.thisN',
         'target_loop.thisIndex', 'target_loop.ran', 'target_loop.order'], axis=1)

# put participant number in personal code
df['Personal Code']=df['Personal Code'].astype('str')
allsubjs=list(df['Personal Code'].unique())
n=1
for i in allsubjs:
    df['Personal Code'][df['Personal Code']==i]=n
    n+=1
del(n,i)

# For every trial, we extract:
trial_RT=[] #Response Times
trial_ID=[] #Trial type (the puppet's ID)
subj=[] #The participant number
X=[] #X location of the mouse click on the screen
Y=[] #Y location of the mouse click on the screen
n=1
for i in range(len(df)-1):
    if i>n:
        if df.iloc[i,4]=='cue_image' and df.iloc[i+1,4]!='cue_image':
        # for every cue, find the next target (expect if we find another cue, we change subject or we hit the end)
            n=i+1
            while df.iloc[n,4]!='target_location' and df.iloc[n,5]==df.iloc[n-1,5] and n<=len(df):
                n+=1
            # if the target was found, add the data relative to the target
            if df.iloc[n,4]=='target_location':
                trial_RT.append(df.iloc[n,3])
                X.append(df.iloc[n,1])
                Y.append(df.iloc[n,2])
                subj.append(df.iloc[n,5])
                # to find the ID, we have to go back until we find green, blue or red
                # except if this is the first trial of a participant, in that case:
                if i==0:
                    trial_ID.append(df.iloc[i,0])
                    first_of_ppt=i
                elif df.iloc[i,5]!=df.iloc[i-1,5]:
                    trial_ID.append(df.iloc[i,0])
                    first_of_ppt=i
                else:
                    m=n-1
                    while df.iloc[m,4]=='target_location' or df.iloc[m,4]=='cue_image':
                        m-=1
                # once the while loop stops, we found a color, which is our ID
                # however, if we went too far back to previous participant, we have to pick the first trial of the participant instead:
                    if df.iloc[i,5]!=df.iloc[m,5]:
                        trial_ID.append(df.iloc[first_of_ppt,0])
                    else:
                        trial_ID.append(df.iloc[m,4]) 
# Trasform the trial ID from colors to numbers
n=-1        
for i in trial_ID:
    n+=1
    if 'blue' in i:
        trial_ID[n]=1
    elif 'red' in i:
        trial_ID[n]=2
    elif 'green' in i:
        trial_ID[n]=3
del(i,n)

# now we have to link the locations reported in X and Y to the valid locations
valid=[[-0.8331310872974386, 0.19987995198079234],
 [-0.8098982070047046, 0.21188475390156064],
 [-0.7773721745948772, 0.22809123649459784],
 [-0.7439168269733403, 0.2424969987995198],
 [-0.7113907945635127, 0.26050420168067223],
 [-0.6742181860951385, 0.2707082833133253],
 [-0.6333283167799268, 0.2755102040816326],
 [-0.5831452953476215, 0.2815126050420168],
 [-0.5431847412441191, 0.28391356542617047],
 [-0.5106587088342917, 0.30012004801920766],
 [-0.4790619916361736, 0.31452581032412963],
 [-0.4335255462624151, 0.3199279711884754],
 [-0.3926356769472033, 0.3247298919567827],
 [-0.3470992315734449, 0.3277310924369748],
 [-0.3043507318348144, 0.33613445378151263],
 [-0.2606729168844746, 0.33853541416566624],
 [-0.2197830475692628, 0.3505402160864346],
 [-0.17424660219550447, 0.35834333733493395],
 [-0.13242741766858335, 0.3565426170468187],
 [-0.0887496027182435, 0.35594237695078035],
 [-0.042283842132775695, 0.3577430972388955],
 [0, 0.36],
 [0.042283842132775695, 0.3577430972388955],
 [0.0887496027182435, 0.35594237695078035],
 [0.13242741766858335, 0.3565426170468187],
 [0.17424660219550447, 0.35834333733493395],
 [0.2197830475692628, 0.3505402160864346],
 [0.2606729168844746, 0.33853541416566624],
 [0.3043507318348144, 0.33613445378151263],
 [0.3470992315734449, 0.3277310924369748],
 [0.3926356769472033, 0.3247298919567827],
 [0.4335255462624151, 0.3199279711884754],
 [0.4790619916361736, 0.31452581032412963],
 [0.5106587088342917, 0.30012004801920766],
 [0.5431847412441191, 0.28391356542617047],
 [0.5831452953476215, 0.2815126050420168],
 [0.6333283167799268, 0.2755102040816326],
 [0.6742181860951385, 0.2707082833133253],
 [0.7113907945635127, 0.26050420168067223],
 [0.7439168269733403, 0.2424969987995198],
 [0.7773721745948772, 0.22809123649459784],
 [0.8098982070047046, 0.21188475390156064],
 [0.8331310872974386, 0.19987995198079234]]

# create function to find the distance between two numbers
def getDist(x1,y1,xy2):
    x2=xy2[0]
    y2=xy2[1]
    dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)  
    return dist  

# for each click, find the valid stimulus location that is closer to the click
location=[]
for i in range(len(X)):
    getmin=[]
    for v in range(len(valid)):
        getmin.append(getDist(X[i],Y[i],valid[v]))
    location.append(getmin.index(min(getmin)))
del(valid,i,v,getmin,m,X,Y,df,first_of_ppt)
# possible additional step: to set as nan the clicks that too far off        
    
# For every trial, we have to also know the actual location of the target
# So we import the actual locations of the targets
stable=pd.read_csv(r'U:\Desktop\S_AdultCuriosity\DataProcessing\stableO.csv',header=None)
variable=pd.read_csv(r'U:\Desktop\S_AdultCuriosity\DataProcessing\volatileO.csv',header=None)
noisy=pd.read_csv(r'U:\Desktop\S_AdultCuriosity\DataProcessing\noisyO.csv',header=None)
stable=stable[0].values.tolist()
noisy=noisy[0].values.tolist()
variable=variable[0].values.tolist()

# in pseudo-randomization 2:
# stable=noisy
# variable=stable
# noisy=variable

# in pseudo-randomization 3:
# stable=variable
# variable=noisy
# noisy=stable

# and we create a column in which we specify the actual location for each target the subjects have seen
actual=[]
count1=-1 #initialize counts
count2=-1
count3=-1
totalcount=-1
todiscard=[] # In case something odd happens, it is reported in this vecor (should be empty)
for i in range(len(trial_ID)): #19 18 18
    totalcount+=1
    if trial_ID[i]==1:
        count1+=1
        if count1<35 and count2<35 and count3<35:
            if subj[i]<22: #Randomization 1
                actual.append(stable[count1])
            if 21<subj[i]<41: #Randomization 2
                actual.append(noisy[count1])
            if subj[i]>40: #Randomization 3
                actual.append(variable[count1])
        else:
            print(totalcount)
            todiscard.append(totalcount)
    elif trial_ID[i]==2:
        count2+=1
        if count1<35 and count2<35 and count3<35:
            if subj[i]<22:
                actual.append(variable[count2])
            if 21<subj[i]<41:
                actual.append(stable[count2])
            if subj[i]>40:
                actual.append(noisy[count2])
        else:
            print(totalcount)
            todiscard.append(totalcount)
    elif trial_ID[i]==3:
        count3+=1
        if count1<35 and count2<35 and count3<35:
            if subj[i]<22:
                actual.append(noisy[count3])
            if 21<subj[i]<41:
                actual.append(variable[count3])
            if subj[i]>40:
                actual.append(stable[count3])
        else:
            print(totalcount)
            todiscard.append(totalcount)
    if i==len(trial_ID)-1:
        break
    if subj[i]!=subj[i+1]:
        count1=-1
        count2=-1
        count3=-1
del(i,stable,noisy,variable,count1,count2,count3)


# now compute the distance between actual and predicted location (i.e, the error)
error=abs(np.array(actual) - np.array(location))

# now compute the learning progress as reduction in prediction error
# so it is the difference between prediction error at time t-1 and time t
# negative values: performance in worsening
# positive values: performance is improving
Progress=["nan"]
for i in range(1,len(error)):
    if subj[i]==subj[i-1]:
        if trial_ID[i]==trial_ID[i-1]:
            Progress.append(error[i-1]-error[i])
        elif trial_ID[i]==trial_ID[i-2]:
            Progress.append(error[i-2]-error[i])
        elif trial_ID[i]==trial_ID[i-3]:
            Progress.append(error[i-3]-error[i])
        else:
            Progress.append("nan")
    else:
        Progress.append("nan")
del(i)

# then we need to know is when participants switched from one character to the other
switch=[]
TimeFromLast=[0]
n=0
for i in range(1,len(trial_ID)):
    if subj[i]==subj[i-1]:
        if trial_ID[i]!=trial_ID[i-1]:
            switch.append(1)
            n=0
        else:
            switch.append(0)
            n+=1
    else:
        switch.append(0)
        n=0
    TimeFromLast.append(n)    
switch.append(0)

#And the novelty of the stimulus
novelty=[]
n1=0
n2=0
n3=0
for i in range(0,len(trial_ID)):
    if subj[i]==subj[i-1]: 
        if trial_ID[i]==1:
            novelty.append(n1)
            n1+=1
        elif trial_ID[i]==2:
            novelty.append(n2)
            n2+=1
        elif trial_ID[i]==3:
            novelty.append(n3)
            n3+=1
    else:
        n1=0
        n2=0
        n3=0
        if trial_ID[i]==1:
            novelty.append(n1)
            n1+=1
        elif trial_ID[i]==2:
            novelty.append(n2)
            n2+=1
        elif trial_ID[i]==3:
            novelty.append(n3)
            n3+=1
# now store everything in one dataframe
df = pd.DataFrame()
df['subj']  = subj
df['char']  = trial_ID
df['RT']  = trial_RT
df['belief']  = location
df['actual']  = actual
df['PE']  = error
df['LP']  = Progress
df['switch'] = switch
df['Time']=TimeFromLast
df['novelty']=novelty
            
# Now we use some of the info preprocessed so far as input for the Nassar et al. (2010) model
# Which runs on matlab. So we export a csv for matlab with xdata, ydata, nn, haz, newBlock
# xdata = outcomes
# ydata = subject estimates
# nn = standard deviation of distribution
nn = [3,2,8] # these were specified when generating the distributions of the three puppets
# haz = hazard rate
haz = [4/35,7/35,3/35] #how many times a puppet switches location, over total number of trials
# newBlock = binary vector, true for each trial that was the first in a new block, false everywhere else
matdf=pd.DataFrame({"belief":[], 
                    "actual":[],
                    "PE":[],
                    "LP":[],
                    "newBlock":[],
                    "nn":[],
                    "haz":[]})
for i in df['subj'].unique():
    thisdf=df[df['subj']==i]
    for n in range(1,4):
        ndf=thisdf[thisdf['char']==n]
        ndf=ndf.drop(['subj','char','RT','switch','Time','novelty'], axis=1)
        newBlock=np.zeros(len(ndf))
        newBlock[0]=1
        thisnn=np.zeros(len(ndf))+nn[n-1]
        thishaz=np.zeros(len(ndf))+haz[n-1]
        ndf['newBlock']=newBlock
        ndf['nn']=thisnn
        ndf['haz']=thishaz
        matdf=matdf.append(ndf)
matdf2=pd.DataFrame()
matdf2["belief"]=matdf["belief"]
matdf2["actual"]=matdf["actual"]
matdf2["newBlock"]=matdf["newBlock"]
matdf2["nn"]=matdf["nn"]
matdf2["haz"]=matdf["haz"]
matdf2.to_csv(r'U:\Desktop\S_AdultCuriosity\BayesianModel\matdata.csv')

########################## HERE MATLAB HAPPENS ###############################