from wand.image import Image
from wand.color import Color
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.units import cm, mm
from reportlab.platypus import Table, TableStyle
from reportlab.lib.pagesizes import A4, cm
from reportlab.lib import colors
from Deminer import *

def DeMinEr2Report(sample,sample_Ctrl1,sample_Ctrl2,gene,region=''):
    ### This function creates a report sample_gene.pdf of mutations called using DeMiner algorithm.
    ### Parameters:    sample: string     
    ###                        File "sample.csv" must be in nucleotideCounts directory
    ###                sample_Ctrl1: string     
    ###                        File "sample_Ctrl1.csv" must be in nucleotideCounts directory  
    ###                sample_Ctrl2: string     
    ###                        File "sample_Ctrl1.csv" must be in nucleotideCounts directory
    ###                gene: string     
    ###                        File "gene.fasta" must be in nucleotideCounts directory
    ###                region: range (optional)
    ###                        Specify the region of gene sequence to be analyzed. Use when sequencing  
    ###                        depth is not uniform on gene (typically 100bp at each end shows lower depth).
    ###                        region must be included in range(length(gene)); 
    ###                        if none is provided, whole gene sequence is analyzed.
    
    if region=='':
        region=len(SeqIO.read(gene+'.fasta','fasta').seq)
    
    # Compute alternative (non-reference) nucleotide frequencies and counts in sample (raw) and both controls
    Test,Testd=getAltFreq(sample,gene)
    CTRL1,CTRL1d=getAltFreq(sample_Ctrl1,gene)
    CTRL2,CTRL2d=getAltFreq(sample_Ctrl2,gene)
    Mutcounts=Test*np.transpose([Testd,Testd,Testd,Testd])

    # Compute processed mutation frequencies using DeMinEr Filter/error-correction steps
    TestF=filterTestvsControl(Test,Testd,CTRL1,CTRL2,sample,gene,region,t=3.29,foldChange=1.5)
    MutcountsF=TestF*np.transpose([Testd,Testd,Testd,Testd])
    
    # Compute mutation profile from processed data
    mutationProfile=computeMutationProfile(MutcountsF,Testd,gene,region)
    
    # Identify AID hotspots on target sequence
    (HotSpot_start_fw,HotSpot_start_rv) = getHotspots(gene)
    HotSpotsregion=[i for i in set(HotSpot_start_fw+HotSpot_start_rv) if i in region]
    
    plotMutations(Test,Testd,gene,region,'Test_'+sample+'_raw')
    plotMutations(TestF,Testd,gene,region,'Test_'+sample+'_filter')
    
    # Report mutations in a pdf file
    pdf = Canvas(sample+'_'+gene+'.pdf', pagesize=A4)
    pdf.setFont('Helvetica', 14)
    pdf.setLineWidth(.3)
    pdf.drawString(3*cm, 27.5*cm, 'DeMinEr analysis of '+gene+' mutations in '+sample)
    pdf.line(2*cm, 27*cm,19*cm, 27*cm,)
    pdf.setFont('Helvetica', 10)
    pdf.setFont('Helvetica-Bold', 10)
        # Mutations from raw data
    pdf.drawString(2*cm, 25.5*cm, '>  Statistics of observed mutations from raw counts')
    pdf.setFont('Helvetica', 10)
    (w,h)=PdftoPng('Test_'+sample+'_raw.pdf')
    pdf.drawImage('Test_'+sample+'_raw.png', 2*cm, 19.5*cm,width=7.5*cm,height=7.5*h/w*cm)
    (w,h)=PdftoPng('Test_'+sample+'_raw_bases.pdf')
    pdf.drawImage('Test_'+sample+'_raw_bases.png', 11*cm, 19.5*cm,width=7.5*cm,height=7.5*h/w*cm)
        # Computation of raw mutation rate and hotspots enrichment
    pdf.drawString(3*cm, 19*cm, '- Mutation rate: '\
                       +str(round(np.sum(Mutcounts[region])/np.sum(np.array(Testd)[region])*1e5,1))+' per 100kb')
    pdf.drawString(3*cm, 18.5*cm, '- Mutation enrichment in AID-hotspots: '\
                       +str(round(np.sum(Mutcounts[HotSpotsregion])/np.sum(Mutcounts[region])/(len(HotSpotsregion)/len(region)),2)))
    pdf.setFont('Helvetica-Bold', 10)
        # Mutations from processed data
    pdf.drawString(2*cm, 17*cm, '>  Statistics of observed mutations from DeMinEr-processed counts')
    pdf.setFont('Helvetica', 10)
    (w,h)=PdftoPng('Test_'+sample+'_filter.pdf')
    pdf.drawImage('Test_'+sample+'_filter.png', 2*cm, 11*cm,width=7.5*cm,height=7.5*h/w*cm)
    (w,h)=PdftoPng('Test_'+sample+'_filter_bases.pdf')
    pdf.drawImage('Test_'+sample+'_filter_bases.png', 11*cm, 11*cm,width=7.5*cm,height=7.5*h/w*cm)
        # Computation of processed mutation rate and hotspots enrichment
    pdf.drawString(3*cm, 10.5*cm, '- Mutation rate: '\
                       +str(round(np.sum(MutcountsF[region])/np.sum(np.array(Testd)[region])*1e5,1))+' per 100kb')
    pdf.drawString(3*cm, 10*cm, '- Mutation enrichment in AID-hotspots: '\
                       +str(round(np.sum(MutcountsF[HotSpotsregion])/np.sum(MutcountsF[region])/(len(HotSpotsregion)/len(region)),2)))
    pdf.setFont('Helvetica-Bold', 10)
        # QC of DeMiner filtering process
    pdf.drawString(2*cm, 8.7*cm, '>  DeMinEr filter')
    pdf.setFont('Helvetica', 10)
    (w,h)=PdftoPng('Filter_'+sample+'.pdf')
    pdf.drawImage('Filter_'+sample+'.png', 2*cm, 3*cm,width=7.5*cm,height=7.5*h/w*cm)
    pdf.drawString(3*cm, 2.5*cm, '- Mean depth (x 1000): '+str(int(round(np.mean(np.array(Testd)[region])/1000,0))))
    pdf.setFont('Helvetica-Bold', 10)
    pdf.drawString(12*cm, 8.7*cm, '>  Mutation profile')
    pdf.setFont('Helvetica', 8)
    t=Table(mutationProfile,5*[0.8*cm], 5*[0.8*cm])
    t.setStyle(TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                            ('FONT',(0,0),(0,-1),'Helvetica-Bold'),
                            ('FONT',(0,0),(-1,0),'Helvetica-Bold'),
                            ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                            ('BOX', (1,0), (-2,-1), 0.25, colors.black),
                            ('BOX', (0,1), (1,-1), 0.25, colors.black),
                               ]))
    width, height = A4
    t.wrapOn(pdf, width, height)
    t.drawOn(pdf, 12*cm,4*cm)
    os.remove('Test_'+sample+'_raw.pdf')
    os.remove('Test_'+sample+'_raw_bases.pdf')
    os.remove('Test_'+sample+'_filter.pdf')
    os.remove('Test_'+sample+'_filter_bases.pdf')
    os.remove('Test_'+sample+'_raw.png')
    os.remove('Test_'+sample+'_raw_bases.png')
    os.remove('Test_'+sample+'_filter.png')
    os.remove('Test_'+sample+'_filter_bases.png')
    os.remove('Filter_'+sample+'.pdf')
    os.remove('Filter_'+sample+'.png')
    pdf.showPage()
    pdf.save() 
    return('Done')



def PdftoPng(pdf):
    ### This function converts pdf images to png files for inclusion in pdf report.
    if os.path.exists(pdf[:-4]+".png"):
        os.remove(pdf[:-4]+".png")
    with Image(filename=pdf, resolution=600) as img:
        (w,h)=img.size
        with Image(width=w, height=h, background=Color("white")) as bg:
            bg.composite(img,0,0)
            bg.save(filename=pdf[:-4]+".png")
    return(w,h)