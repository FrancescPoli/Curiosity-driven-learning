# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:32:15 2020

@author: Francesco
"""
# Run this script after S_AdultCuriosity_Preprocessing.py and loopover.m

#Import alpha as computed by the Bayesian model to estimate expected LP
pexp=[]
expLP=[]
n1=0
n2=0
n3=0
alphaavg=pd.read_csv(r'U:\Desktop\S_AdultCuriosity\BayesianModel\alphaavg.csv',header=None)
alphaavg=list(alphaavg.iloc[:,0])

groupalpha=[]
for i in range(0,len(df["subj"])):
    if df["char"][i]==1:
        groupalpha.append(alphaavg[df["Time"][i]])
    elif df["char"][i]==2:
        groupalpha.append(alphaavg[df["Time"][i]+35])
    elif df["char"][i]==3:
        groupalpha.append(alphaavg[df["Time"][i]+70])
        
for i in range(0,len(trial_ID)):
    if subj[i]==subj[i-1]: 
        if trial_ID[i]==1:
            pexp1=pexp1+groupalpha[i]*(error[i]-pexp1)
            pexp.append(pexp1)
            expLP.append(error[i]-pexp1)
            n1+=1
        elif trial_ID[i]==2:
            pexp2=pexp2+groupalpha[i]*(error[i]-pexp2)
            pexp.append(pexp2)
            expLP.append(error[i]-pexp2)
            n2+=1
        elif trial_ID[i]==3:
            pexp3=pexp3+groupalpha[i]*(error[i]-pexp3)
            pexp.append(pexp3)
            expLP.append(error[i]-pexp3)
            n3+=1
    else:
        n1=0
        n2=0
        n3=0
        pexp1=0
        pexp2=0
        pexp3=0
        expLP.append('nan')
        if trial_ID[i]==1:
            pexp.append(pexp1)
            n1+=1
        elif trial_ID[i]==2:
            pexp.append(pexp2)
            n2+=1
        elif trial_ID[i]==3:
            pexp.append(pexp3)
            n3+=1


df['alpha']=groupalpha
df['expLP']=expLP
df['expPE']=pexp
del(Progress,TimeFromLast,actual,alphaavg,error,expLP,groupalpha,haz,i,location,matdf,matdf2,n,n1,n2,n3,ndf,newBlock,nn,novelty,pexp,pexp1,pexp2,pexp3,subj,switch,thisdf,thishaz,thisnn,todiscard,totalcount,trial_ID,trial_RT)

# now structure data for analysis of delta values between puppets
# when a switch is made, we have to know the difference in KL and novelty between the two new characters

deltaEIG=[] # difference in LP
deltaNov=[] # difference in novelty
chosenoption=[]
subj=[]
for s in df['subj'].unique():
    EIG=[-1,-1,-1]
    Nov=[-1,-1,-1]
    thisdf=df[df['subj']==s]
    thisdf=thisdf.reset_index()
    thisdf['expLP'] = thisdf['expLP'].replace(np.inf, "nan")
    thischar=thisdf['char'][0]
    for i in range(len(thisdf['expLP'])): #+thisdf.index.values.astype(int)[0]
        #when there is a change in character
        if thisdf['char'][i] != thischar:
            oldchar=thischar
            thischar=thisdf['char'][i]            
            n=i
            stop=0
            while stop==0 and thisdf['expLP'][n] == "nan":
                n-=1
                if n<0:
                    stop=1
                    n=0
            if n>0:
                # we should check we didn't go too far backwards, to another chartype
                #set the EIG of thischar
                EIG[oldchar-1]=thisdf['expLP'][n]
                Nov[oldchar-1]=thisdf['novelty'][n]
                if EIG[0]!=-1 and EIG[1]!=-1 and EIG[2]!=-1: # when we have all the EIGs
                # we compute the difference in eigs:
                #fisrt we delete the eig of the current char
                    thisEIG=[EIG[0],EIG[1],EIG[2]]
                    thisEIG.pop(oldchar-1)
                    thisNov=[Nov[0],Nov[1],Nov[2]]
                    thisNov.pop(oldchar-1)
                    # if they don't have the same eig
                    if thisEIG[0]!=thisEIG[1]:
                        #high=np.argmax(thisEIG)
                        # now we determine what was option 1 and whether subjects chose it
                        if oldchar == 1:
                            option1=2
                        else:
                            option1=1
                        # now we can check whether they chose option 1
                        if option1 == thisdf['char'][i]:
                            chose1=0
                        else:
                            chose1=1
                        #now we set the delta EIG and the chosen option
                        deltaEIG.append(thisEIG[1]-thisEIG[0])
                        deltaNov.append(-1*(thisNov[1]-thisNov[0])) #multiplied by -1 because high scores = low novelty
                        chosenoption.append(chose1)
                        subj.append(s)

#store everything in new dataframe
EIGdat=pd.DataFrame()                       
EIGdat['deltaEIG']=deltaEIG
EIGdat['choice']=chosenoption
EIGdat['subject']=subj
EIGdat['DeltaNovelty']=deltaNov

del(Nov,chosenoption,deltaEIG,deltaNov,s,subj,thisdf,chose1,i,n,oldchar,option1,stop,thisEIG,thisNov,thischar)

# Now let's deal with the questionnaires
poll=pd.read_csv(r'U:\Desktop\S_AdultCuriosity\DataProcessing\data\results-survey.csv',delimiter=';')
# Now we can get the descriptives for participants
descriptives=pd.DataFrame()
descriptives['subjmean']=np.mean(poll.iloc[:,5])
descriptives['subjstd']=np.std(poll.iloc[:,5])
descriptives['F']=sum(poll.iloc[:,6]=='Female')
descriptives['tot']=len(poll)

# now remove subjects that did something wrong
subj_toexclude=[]
number_toexclude=[]
for i in range(len(poll)):
    if poll.iloc[i,4]=='Yes' or poll.iloc[i,3]=='Randomly':
        subj_toexclude.append(poll.iloc[i,0])     
        
# If participants are in subj_toexclude, they must be excluded from the analysis.
for i in poll.iloc[:,0]:
    if i in subj_toexclude:
        number_toexclude.append(poll[poll['Subject']==i].index.values.astype(int)[0]+1)
        poll=poll.drop(poll.index[poll.iloc[:,0]==i].tolist(),axis=0)
del(subj_toexclude)

# also from the df and EIGdata analyses!
for i in number_toexclude:
    df=df.drop(index=df[df['subj'] == i].index)
del(allsubjs,i,number_toexclude)

# Now we can get the descriptives for participants aftre exclusion
descriptives=pd.DataFrame()
descriptives['subjmean']=np.mean(poll.iloc[:,5])
descriptives['subjstd']=np.std(poll.iloc[:,5])
descriptives['F']=sum(poll.iloc[:,6]=='Female')
descriptives['tot']=len(poll)

# get favourite and least favourite 
polldata = pd.DataFrame()
count=-1
chartype_most=[]
chartype_least=[]
for i in range(1,len(poll['Subject'])+1):
    count+=1
    if poll.iloc[count,1]=='Blue':
        if i < 22:
            chartype_most.append(1)
        elif 21 < i < 41:
            chartype_most.append(3)
        elif i > 40:
            chartype_most.append(2)
    elif poll.iloc[count,1]=='Red':
        if i < 22:
            chartype_most.append(2)
        elif 21 < i < 41:
            chartype_most.append(1)
        elif i > 40:
            chartype_most.append(3)
    elif poll.iloc[count,1]=='Green':
        if i < 22:
            chartype_most.append(3)
        elif 21 < i < 41:
            chartype_most.append(2)
        elif i > 40:
            chartype_most.append(1)
    else:
        chartype_most.append("nan")

count=-1
for i in range(1,len(poll['Subject'])+1):
    count+=1
    if poll.iloc[count,2]=='Blue':
        if i < 22:
            chartype_least.append(1)
        elif 21 < i < 41:
            chartype_least.append(3)
        elif i > 40:
            chartype_least.append(2)
    elif poll.iloc[count,2]=='Red':
        if i < 22:
            chartype_least.append(2)
        elif 21 < i < 41:
            chartype_least.append(1)
        elif i > 40:
            chartype_least.append(3)
    elif poll.iloc[count,2]=='Green':
        if i < 22:
            chartype_least.append(3)
        elif 21 < i < 41:
            chartype_least.append(2)
        elif i > 40:
            chartype_least.append(1)
    else:
        chartype_least.append("nan")

polldata['chartype_most']=chartype_most
polldata['chartype_least']=chartype_least  

# check mean and std of trail number in participants
np.mean(df['subj'].value_counts())
np.std(df['subj'].value_counts())

#export all data
# export as csv for statistical analysis on R both df as data and EIGdat as EIGdata
df.to_csv(r'U:\Desktop\S_AdultCuriosity\DataAnalysis\datalp.csv')
EIGdat.to_csv(r'U:\Desktop\S_AdultCuriosity\DataAnalysis\EIGdata.csv')
polldata.to_csv(r'U:\Desktop\S_AdultCuriosity\DataAnalysis\polldata.csv')

# PLOTS 

# PLOT FOR STIMULI DISTRIBUTION    

# We import the actual locations of the targets
stable=pd.read_csv(r'C:\Users\Francesco\Desktop\Hide&Seek_Analysis\stableO.csv',header=None)
variable=pd.read_csv(r'C:\Users\Francesco\Desktop\Hide&Seek_Analysis\volatileO.csv',header=None)
noisy=pd.read_csv(r'C:\Users\Francesco\Desktop\Hide&Seek_Analysis\noisyO.csv',header=None)
stable=stable[0].values.tolist()
noisy=noisy[0].values.tolist()
variable=variable[0].values.tolist()

idealmat=pd.read_csv(r'C:\Users\Francesco\Desktop\Hide&Seek_Analysis\idealmat.csv',header=None)
idealmat.columns= ['volatileM','volatileO','Bvolatile','Sigvolatile','Alphavolatile', 'stableM','stableO','Bstable','Sigstable','Alphastable', 'noisyM','noisyO','Bnoisy','Signoisy','Alphanoisy'];

os.chdir(r"C:\Users\Francesco\Desktop\Hide&Seek_Analysis")
plt.figure(figsize=(12, 4)) 
 
plt.subplot(131)
plt.plot(idealmat['Bstable'],c='grey')
plt.plot(idealmat['stableO'],c='#4c72b0')
plt.ylim(-1, 44)
plt.ylabel('Location')
plt.xlabel('Trial')
ax = plt.gca()
ax.title.set_text('Intermediate Noise\nIntermediate Volatility')
ax.set_facecolor('xkcd:white')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

plt.subplot(132)
plt.plot(idealmat['Bvolatile'],c='grey')     
plt.plot(idealmat['volatileO'],c='#c44e52')
plt.ylim(-1, 44)
plt.xlabel('Trial')
ax = plt.gca()
ax.title.set_text('Low Noise\nHigh Volatility')
ax.set_facecolor('xkcd:white')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

plt.subplot(133)
plt.plot(idealmat['Bnoisy'],c='grey')
plt.plot(idealmat['noisyO'], c='#55a868')
plt.ylim(-1, 44)
plt.xlabel('Trial')
ax = plt.gca()
ax.title.set_text('High Noise\nLow Volatility')
ax.set_facecolor('xkcd:white')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
#plt.show()        

plt.savefig('chars.png',dpi=300)
 

# Plot frequencies Explicit answers
data = [[np.sum(polldata['chartype_most']==1), np.sum(polldata['chartype_least']==1)],
[np.sum(polldata['chartype_most']==2), np.sum(polldata['chartype_least']==2)],
[np.sum(polldata['chartype_most']==3), np.sum(polldata['chartype_least']==3)]]
X = np.arange(2)
fig = plt.figure() 
ax = fig.add_axes([0,0,1,1])
ax.set_facecolor('xkcd:white')
ax.bar(X + 0.00, data[0], color = 'b', width = 0.25)
ax.bar(X + 0.25, data[1], color = 'r', width = 0.25)
ax.bar(X + 0.50, data[2], color = 'g', width = 0.25)
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
plt.ylabel('Frequency')
plt.legend(['Intermediate', 'High Volatility', 'High Noise'])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[2] = 'Most Favourite'
labels[6] = 'Least Favourite'
ax.set_xticklabels(labels)

plt.yticks(np.arange(0, 19, 5))

# plot AIC and BIC for the three models + null model
data=[[3008, 3550, 3572, 3575],
[3306, 3857, 3900, 3877]]
X = np.arange(4)
fig = plt.figure(figsize=(6, 4))
ax = fig.add_axes([0,0,1,1])
ax.set_facecolor('xkcd:white')
ax.bar(X + 0.00, data[0], color = '#eb6600', width = 0.4)
ax.bar(X + 0.4, data[1], color = '#008db0', width = 0.4)
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.set_ylim((2500, 4000)) 
plt.legend(['AIC', 'BIC'])
#plt.ylabel('Frequency')

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[1] = '             Learning Progress'
labels[3] = '               Novelty'
labels[5] = '               Random Search'
labels[7] = '               Null Model'
ax.set_xticklabels(labels)

plt.yticks(np.arange(2500, 4001, 500))
plt.savefig('AIC.png',dpi=300)


data=[[2934],
[3478],
[3487],
[3489],]
X = np.arange(1)
fig = plt.figure(figsize=(6, 4))
ax = fig.add_axes([0,0,1,1])
ax.set_facecolor('xkcd:white')
ax.bar(X + 0.00, data[0], color = '#d1e6fa', width = 0.4)
ax.bar(X + 0.5, data[1], color = '#d1e6fa', width = 0.4)
ax.bar(X + 1, data[2], color = '#d1e6fa', width = 0.4)
ax.bar(X + 1.5, data[3], color = '#d1e6fa', width = 0.4)
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')
ax.set_ylim((2600, 3600)) 
#plt.legend(['AIC', 'BIC'])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[2] = 'Learning Progress'
labels[4] = 'Novelty'
labels[6] = 'Null Model'
labels[8] = 'Random Search'
ax.set_xticklabels(labels)
plt.ylabel('AIC', fontsize=18)
plt.xticks(fontsize=13)

plt.yticks(np.arange(2600, 3601, 200))
plt.savefig('AIC.png',dpi=300)