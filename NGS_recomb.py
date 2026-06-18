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
recomb_path=home_path+"recomb_files/"

## script gets samplenames by checking whats in the sample list

file=open(home_path+"samples.csv","r")
csv_file=csv.reader(file, delimiter=";")

samplelist=[]
samplename=""
reads=""

depth_av=[]

trimmoutput=""
trimmomatic_options="ILLUMINACLIP:/home/hoeflet95/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

for line in csv_file:                                               #loading samplenames into samplename list
    name=str(line[2])+"_"+str(line[5])+"_P"+str(line[3])+"Rep"+str(line[4])
    samplelist.append(line)
    
    ## trimming raw reads to get rid of leftover adapter sequences with Trimmomatic    
        
    reads=raw_path+line[0]+" "+raw_path+line[1]+" "
    trimmoutput=trimmed_path+name+"_forward_paired.fq.gz "+trimmed_path+name+"_forward_unpaired.fq.gz "+trimmed_path+name+"_reverse_paired.fq.gz "+trimmed_path+name+"_reverse_unpaired.fq.gz "
    os.system("java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE "+reads+trimmoutput+trimmomatic_options)

## from here on we work only with our trimmed samplenames, the ones without Illumina samplenumbers

file.close()
os.system("rm "+trimmed_path+"*unpaired*") #removing unpaired reads

       

for line in samplelist:  
    name=str(line[2])+"_"+str(line[5])+"_P"+str(line[3])+"Rep"+str(line[4])
    
    ## Starting recombination analysis with PEAR and ViReMa
    os.system("pear -f "+trimmed_path+name+"_forward_paired.fq.gz -r "+trimmed_path+name+"_reverse_paired.fq.gz -o "+trimmed_path+name)
 
    
    os.system("python "+home_path+"ViReMa.py "+reference_path+"HSV-1_F_trimmed.fa "+trimmed_path+name+".assembled.fastq "+name+".sam -BED --Seed 20 --Aligner bwa --p 15 --Output_Tag "+name+" --Output_Dir "+recomb_path+name)
    

print("done!")
    
    
    
    
