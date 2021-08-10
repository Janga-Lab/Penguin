
#note user should provide the coordinate file and the fastq file
import subprocess
import sys
import argparse

#https://stackoverflow.com/questions/21579892/using-argparse-to-for-filename-input 
 
# Initialize parser
parser = argparse.ArgumentParser("start of getting RNA fastq files")

#reference file argument
parser.add_argument('-r', action='store', dest='ref_input', help='Provide reference genome file')

#fastq file argument 
parser.add_argument('-f', action='store', dest='fastq_input', help='Provide Fastq file')
 

# Read arguments from command line
args = parser.parse_args()

ref = args.ref_input
fastq_reads = args.fastq_input


##################################################################################
#a=sys.argv[1]
#b=sys.argv[2]
f = open("hek.sam", "w")
f1 = open("hek.bam", "w")
f3 = open("hek-reads-ref.eventalign.txt", "w")
#f4 = open("Pseu_Modification_coors_ns_hek.txt", "w") #coordinate file
####subprocess.run(["minimap2", "-ax", "map-ont", "--split-prefix", "/tmp/temp_name", "ref.fa", "reads.fastq"], stdout=f)    #minimap2 -ax map-ont --split-prefix /tmp/temp_name  ref.fa  reads.fastq > hek.sam
subprocess.run(["minimap2", "-ax", "map-ont", "--split-prefix", "/tmp/temp_name", ref, fastq_reads], stdout=f) 
subprocess.run(["nanopolish" ,"index", "-d", "fast5_files/", "reads.fastq"])   #nanopolish index -d fast5_files/ reads.fastq 
subprocess.run(["samtools" ,"view", "-S", "-b", "hek.sam"],stdout=f1)#samtools view -S -b hela.sam > hela.bam 
subprocess.run(["samtools" ,"sort", "hek.bam", "-o", "hek.sorted.bam"]) #samtools sort hela.bam -o hela.sorted.bam
subprocess.run(["samtools" ,"index", "hek.sorted.bam"]) #samtools index hela.sorted.bam
subprocess.run(["samtools","quickcheck" ,"hek.sorted.bam"]) #$samtools quickcheck hela.sorted.bam
subprocess.run(["nanopolish","eventalign" ,"--reads", "reads.fastq", "--bam", "hek.sorted.bam","--genome", "ref.fa", "--scale-events"], stdout=f3) #$nanopolish eventalign  --reads reads.fastq --bam  hela.sorted.bam  --genome ref.fa --scale-events > hela-reads-ref.eventalign.txt
subprocess.run(["python","flash_m5c_coors.py"])
subprocess.run(["python","SVM_onehot.py",])





