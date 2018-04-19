import os
import re
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def getAltFreq(sample,gene):
    ### This function returns NumPy arrays of alternative (non-reference) nucleotide frequencies, 
    ### depth and alternative nucleotide counts.
    ### Parameters:    sample: string     
    ###                        File "sample.csv" must be in nucleotideCounts directory
    ###                gene: string     
    ###                        File "gene.fasta" must be in nucleotideCounts directory    
    ### Outputs:       Freqs: Numpy array [length(gene) x 4]
    ###                        Alternative nucleotide frequencies at each position
    ###                Depths: list  [length(gene)]   
    ###                        Depth at each position    
    ###                Mutcounts: Numpy array   [length(gene) x 4]   
    ###                        Alternative nucleotide counts at each position      
    
    # Load nucleotide counts csv file
    csv='nucleotideCounts/'+sample+'_'+gene+'.csv'
    Counts=[]
    with open(csv,'r') as f:
        for line in f:
            Counts.append([int(c) for c in line.rstrip(';\n').split(';')])
    
    # Compute nucleotide frequencies and depth
    Depths=[]
    Freqs=[]
    for count in Counts:
        Depths.append(sum(count))
        Freqs.append([c/sum(count) for c in count])
    Counts=np.array(Counts)
    Freqs=np.array(Freqs)
    
    # Filter reference nucleotide at each position
    (n,m)=np.shape(Freqs)
    for i in range(n):
        for j in range(m):
            if Freqs[i,j]>0.5:
                Freqs[i,j]=0
    return(Freqs,Depths)

def filterTestvsControl(Test,Testd,CTRL1,CTRL2,sample,gene,region,t=3.29,foldChange=1.5):
    ### This function returns NumPy arrays of prediction intervals.
    ### Parameters:    Test: Numpy array [length(gene) x 4]     
    ###                        Alternative nucleotide frequencies in sample as returned by getAltFreq
    ###                Testd: Numpy array [length(gene) x 1]     
    ###                        Position-wise depth in sample as returned by getAltFreq
    ###                CTRL1: Numpy array [length(gene) x 4]     
    ###                        Alternative nucleotide frequencies in Ctrl1 as returned by getAltFreq
    ###                CTRL2: Numpy array [length(gene) x 4]     
    ###                        Alternative nucleotide frequencies in Ctrl2 as returned by getAltFreq
    ###                sample: string
    ###                gene: string     
    ###                        File "gene.fasta" must be in nucleotideCounts directory   
    ###                region: range
    ###                        Specify the region of gene sequence.
    ###                t: float (optional)
    ###                        Student's distribution coefficient used in PI computation
    ###                foldChange: float (optional)
    ###                        foldChange threshold in Test vs. Control comparison
    
    ### Outputs:       TestF: Numpy array [length(gene) x 4]
    ###                        DeMinEr-corrected alternative nucleotide frequencies in test sample
    
    
    # Compute prediction intervals of Ctrl1 vs. Ctrl2 correlation
    d_PIF=np.zeros_like(CTRL1)
    x_PI=np.zeros([100,4])
    d_PI=np.zeros([100,4])
    for i in range(4):
        n=len(region)
        m=np.mean(CTRL1[region,i])
        s=sum(np.power(CTRL2[region,i]-CTRL1[region,i],2))
        mse=np.sqrt(sum(np.power(CTRL2[region,i]-CTRL1[region,i],2))/(n-2))
        x_PI[:,i]=np.arange(0,0.1,0.001)
        d_PI[:,i]=t*mse*np.sqrt(1+1/n+np.power(x_PI[:,i]-m,2)/s)
        d_PIF[:,i]=t*mse*np.sqrt(1+1/n+np.power(CTRL1[:,i]-m,2)/s)

    # Filter off-correlation mutations 
    Filter=Test>CTRL1+d_PIF
    FilterFC=Test>foldChange*CTRL1
    TestF=(Test-CTRL1)*Filter*FilterFC
    
    # Plot 'Test vs. Ctrl1 nucleotide frequencies / highlight filtered events
    colors=['#57bddb','#ea4848','#3bb971','#f39c12']
    for i in range(4):
        plt.scatter(CTRL1[region,i],Test[region,i],color=colors[i],alpha=0.2)
        plt.plot([0,1],[0,1],'g')
        plt.plot(x_PI[:,i],x_PI[:,i]+d_PI[:,i],color=colors[i])
        plt.plot(x_PI[:,i],x_PI[:,i]-d_PI[:,i],color=colors[i])
        plt.scatter(CTRL1[region,i]*Filter[region,i]*FilterFC[region,i],
                    Test[region,i]*Filter[region,i]*FilterFC[region,i],color='r')    
        plt.rcParams.update({'font.size': 14})
        plt.axis([0,0.01,0,0.01])
        plt.xlabel('Mutation frequency in Ctrl')
        plt.ylabel('Mutation frequency in Test')
    plt.savefig('Filter_'+sample+'.pdf',bbox_inches='tight')
    plt.close()
    return(TestF)

def plotMutations(Test,Testd,gene,region,figName,display=False):
    ### This function returns NumPy arrays of prediction intervals.
    ### Parameters:    Test: Numpy array [length(gene) x 4]     
    ###                        Alternative nucleotide frequencies in sample as returned by getAltFreq
    ###                Testd: Numpy array [length(gene) x 1]     
    ###                        Position-wise depth in sample as returned by getAltFreq
    ###                gene: string     
    ###                        File "gene.fasta" must be in nucleotideCounts directory   
    ###                region: range
    ###                        Specify the region of gene sequence.
    ###                figName: string     
    ###                        File name to save generated figures   
    
    
    Mutcounts=Test*np.transpose([Testd,Testd,Testd,Testd])
    colors=['#57bddb','#ea4848','#3bb971','#f39c12']
    # Find hotspots' position
    refSeq=str(SeqIO.read(gene+'.fasta','fasta').seq)
    (HotSpot_start_fw,HotSpot_start_rv) = getHotspots(gene)
    
    # Plot position-wise mismatch frequency
    plt.plot(Test[:,0]+Test[:,1]+Test[:,2]+Test[:,3],'b')
    plt.rcParams.update({'font.size': 14})
    plt.axis([0,len(refSeq),0,0.02])
    plt.xlabel('Position (bp)')
    plt.ylabel('Mutation frequency')
    ax=plt.gca()
    for i in range(0,len(HotSpot_start_fw),4):
        ax.add_patch(patches.Rectangle((HotSpot_start_fw[i], 0), 4, 1,color='m',alpha=0.1))
    for i in range(0,len(HotSpot_start_rv),4):
        ax.add_patch(patches.Rectangle((HotSpot_start_rv[i], 0), 4, 1,color='g',alpha=0.1))
    if display:
        plt.show()
    plt.savefig(figName+'.pdf',bbox_inches='tight')
    plt.close()
    
    # Plot position-wise mismatch frequency for each nucleotide
    ymax=0.02
    f, axarr = plt.subplots(2, 2)
    plt.rcParams.update({'font.size': 10})
    axarr[0, 0].plot(Test[:,0],color=colors[0])
    axarr[0, 0].set_ylim([0,ymax])
    axarr[0, 0].set_xlim([0,len(refSeq)])
    axarr[0, 0].set_title('Mut. to A = '\
                          +str(round(np.sum(Mutcounts[region,0])/np.sum(np.array(Testd)[region])*1e5,2))+'/100kb',
                         color=colors[0])
    axarr[0, 1].plot(Test[:,1],color=colors[1])
    axarr[0, 1].set_ylim([0,ymax])
    axarr[0, 1].set_xlim([0,len(refSeq)])
    axarr[0, 1].set_title('Mut. to C = '\
                          +str(round(np.sum(Mutcounts[region,1])/np.sum(np.array(Testd)[region])*1e5,2))+'/100kb',
                         color=colors[1])
    axarr[1, 0].plot(Test[:,2],color=colors[2])
    axarr[1, 0].set_ylim([0,ymax])
    axarr[1, 0].set_xlim([0,len(refSeq)])
    axarr[1, 0].set_title('Mut. to G = '\
                          +str(round(np.sum(Mutcounts[region,2])/np.sum(np.array(Testd)[region])*1e5,2))+'/100kb',
                         color=colors[2])
    axarr[1, 1].plot(Test[:,3],color=colors[3])
    axarr[1, 1].set_ylim([0,ymax])
    axarr[1, 1].set_xlim([0,len(refSeq)])
    axarr[1, 1].set_title('Mut. to T = '\
                          +str(round(np.sum(Mutcounts[region,3])/np.sum(np.array(Testd)[region])*1e5,2))+'/100kb',
                         color=colors[3])
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    axarr[1, 1].set_xlabel('Position (bp)')
    axarr[1, 0].set_xlabel('Position (bp)')
    axarr[1, 0].set_ylabel('Mut. frequency')
    axarr[0, 0].set_ylabel('Mut. frequency')
    if display:
        plt.show()
    plt.savefig(figName+'_bases.pdf',bbox_inches='tight')
    plt.close()
    return


def computeMutationProfile(MutcountsF,Testd,gene,region):
    ### This function returns NumPy arrays of prediction intervals.
    ### Parameters:    MutCountsF: Numpy array [length(gene) x 4]     
    ###                        Peocessed alternative nucleotide frequencies in sample
    ###                Testd: Numpy array [length(gene) x 1]     
    ###                        Position-wise depth in sample as returned by getAltFreq
    ###                gene: string     
    ###                        File "gene.fasta" must be in nucleotideCounts directory   
    ###                region: range
    ###                        Specify the region of gene sequence.    
    ### Outputs:       mutationProfile: array
    ###                        Composition-corrected conversion table
    refSeq=str(SeqIO.read(gene+'.fasta','fasta').seq)
    alphabet={'A':0,'C':1,'G':2,'T':3}
    alphabet2={'A':0,'T':1,'C':2,'G':3}
    
    # Build conversion table
    mutations=np.zeros([4,5])
    for i in region:
        for j in range(4):
            mutations[alphabet[refSeq[i]],j]+=MutcountsF[i,j]
    temp=np.zeros_like(mutations)
    temp[:,0]=mutations[[0,3,1,2],0]
    temp[:,1]=mutations[[0,3,1,2],3]
    temp[:,2]=mutations[[0,3,1,2],1]
    temp[:,3]=mutations[[0,3,1,2],2]
    mutations=temp
    
    # Apply sequence composition-correction and normalize
    sumBases=np.zeros(4)
    for i in region:
        sumBases[alphabet2[refSeq[i]]]+=Testd[i]
    corr=0.25/(sumBases/np.sum(sumBases))
    mutationsCorr=np.zeros_like(mutations)
    for i in range(4):
        for j in range(5):
            mutationsCorr[i,j]=mutations[i,j]*corr[i]
    S=max(1,np.sum(mutationsCorr))
    for i in range(4):
        for j in range(5):
            mutationsCorr[i,j]=round(mutationsCorr[i,j]/S*100)
    for i in range(4):
        mutationsCorr[i,4]=sum(mutationsCorr[i,:])
    mutationProfile=[[' ','A','T','C','G',' ']]
    for i in range(4):
        mutationProfile.append(['ATCG'[i]]+[int(k) for k in mutationsCorr[i,:]])
    return(mutationProfile)


def getHotspots(gene):
    ### This function returns hotspots' position in gene sequence.
    refSeq=str(SeqIO.read(gene+'.fasta','fasta').seq)
    HotSpot_start_fw=[]
    HotSpot_start_rv=[]
    for m in re.finditer("[AT][AG]C[CT]", refSeq):
        HotSpot_start_fw.append(m.start()+2)
    for m in re.finditer("[AG]G[CT][AT]", refSeq):
        HotSpot_start_rv.append(m.start()+1)
    return(HotSpot_start_fw,HotSpot_start_rv)
