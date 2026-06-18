import os 
import csv

## defining sample and data paths, important are:
## home_path, raw_path (raw reads are located there), timmed_path (that's where trimmed reads are saved),
## reference_path (indexed references are saved there), vcf_path (variants calls are saved there),
## bam_path (where mapping alignments are saved), coverage_path (where coverage data is saved) and 
## recomb_path (where ViReMa output is saved)

home_path="/path/to/home_directory/"

raw_path=home_path+"raw_data/"
trimmed_path=home_path+"trimmed_reads/"
reference_path=home_path+"reference/"
vcf_path=home_path+"vcf_files/"
bam_path=home_path+"bam_files/"
coverage_path=home_path+"coverage/"
recomb_path=home_path+"recomb_files/"

## script gets samplenames by checking whats in the raw read directory

os.system("ls "+raw_path+" > "+home_path+"temp.txt") #getting samplenames
file=open(home_path+"temp.txt","r")

samplelist=[]
samplename=""
samplelist2=[]
reads=""

depth_EGFP=[]
depth_av=[]

trimmoutput=""
trimmomatic_options="ILLUMINACLIP:/home/hoeflet95/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

for line in file:                                               #loading samplenames into samplename list
    samplename=line[:line.find("_R")]
    if samplename not in samplelist:
        samplelist.append(samplename)

file.close()
os.system("rm "+home_path+"temp.txt")

print("trimming samplenames ...")

for name in samplelist:                                         
    samplename=name[:name.find("_")]                             #trimming samplenames
    if samplename not in samplelist2:
        samplelist2.append(samplename)

    ## trimming raw reads to get rid of leftover adapter sequences with Trimmomatic    
        
    reads=raw_path+name+"_R1_001.fastq.gz "+raw_path+name+"_R2_001.fastq.gz "
    trimmoutput=trimmed_path+samplename+"_forward_paired.fq.gz "+trimmed_path+samplename+"_forward_unpaired.fq.gz "+trimmed_path+samplename+"_reverse_paired.fq.gz "+trimmed_path+samplename+"_reverse_unpaired.fq.gz "
    os.system("java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE "+reads+trimmoutput+trimmomatic_options)

## from here on we work only with our trimmed samplenames, the ones without Illumina samplenumbers

os.system("rm "+trimmed_path+"*unpaired*") #removing unpaired reads     

for name in samplelist2:  
    
    os.system("cd "+home_path)
    
    ## mapping trimmed paired reads to reference with Burrows-Wheeler Aligner, reference is already indexed!
    ## SAM file is sorted as well as converted to BAM file and indexed afterwards
    
    os.system("bwa mem -t 16 "+reference_path+"HSV-1_F_trimmed.fa "+trimmed_path+name+"_forward_paired.fq.gz "+trimmed_path+name+"_reverse_paired.fq.gz > "+bam_path+name+".sam")  
    os.system("samtools view -Sb "+bam_path+name+".sam | samtools sort - > "+bam_path+name+".bam")
    os.system("rm "+bam_path+name+".sam")
    os.system("samtools index "+bam_path+name+".bam")

    ## Variant call is performed with LoFreq and BCFtools
    ## LoFreq indelqual adds indel quality values to the bam alignment, imporant for indel calls
    
    os.system("lofreq indelqual --dindel -f "+reference_path+"HSV-1_F_trimmed.fa "+bam_path+name+".bam > "+bam_path+name+"_indel.bam")
    os.system("samtools index "+bam_path+name+"_indel.bam")
    os.system("lofreq call-parallel --pp-threads 16 --call-indels -f "+reference_path+"HSV-1_F_trimmed.fa -o "+vcf_path+name+"_lofreq.vcf "+bam_path+name+"_indel.bam")  
    os.system("bgzip "+vcf_path+name+"_lofreq.vcf")
    os.system("tabix "+vcf_path+name+"_lofreq.vcf.gz")

    os.system("bcftools mpileup -f "+reference_path+"HSV-1_F_trimmed.fa "+bam_path+name+"_indel.bam | bcftools call -cv -Ov --ploidy 1 -o "+vcf_path+name+"_c.vcf")
    os.system("bgzip "+vcf_path+name+"_c.vcf")
    os.system("tabix "+vcf_path+name+"_c.vcf.gz")
    
    ## Getting coverage data with SAMtools depth
    
    os.system("samtools depth -r HSV-1_F_trimmed "+bam_path+name+"_indel.bam > "+coverage_path+name+"_cov_av.cov")
    
    ## Converting bam to cram, saving disc space
    
    os.system("samtools view -T "+reference_path+"HSV-1_F_trimmed.fa -C -o "+bam_path+name+".cram "+bam_path+name+"_indel.bam")
    os.system("samtools index "+bam_path+name+".cram")
    os.system("rm "+bam_path+"*.bam*")
    
    ## Calculating average coverage
    
    file_av=open(coverage_path+name+"_cov_av.cov","r")
    coverage_av=[]
        
    for av in file_av:
        coverage_av.append(av.split())
        
    cov_av=0.0

    for av in coverage_av:
        cov_av=cov_av+float(av[2])
        
    if len(coverage_av)>0:    
        cov_av=cov_av/len(coverage_av)
    else:
        cov_av=0    
    
        
  
    depth_av.append([name,cov_av])
    file_av.close()
    
## Writing the average coverage for all samples onto a CSV file
    
out_av=["name","average depth"]

outfile_av=open(coverage_path+"depth_av.csv","w")
out_csv=csv.writer(outfile_av,delimiter=";")
out_csv.writerow(out_av)
for av in depth_av:
    out_csv.writerow(av)

outfile_av.close()

    
print("done!")
    
    
    
    
