#!/usr/bin/bash

####################
#Installations
conda install entrez-direct
conda install sra-tool kit (#or from source code)
conda install -c bioconda parallel-fastq-dump
####################

#SRA download
#query ID
sra_id=SRX7805363 #single end
sra_id=SRR17062757 #paired end
sra_tools=sratoolkit.3.0.1-mac64/bin

#query SRX ID to generate SRR ID for fetching SRA files
esearch -db sra -query SRX7805363 | efetch -format runinfo| cut -d ',' -f 1 | grep SRR

#prefetch saves time and fastens the process
time prefetch SRR11184871
time ${sra_tools}/prefetch  SRR17062757
prefetch $(<SRA_Acc_List.txt) #list of IDs

#check file integrity and check md5sums
vdb-validate SRR11184871

#fasterq with multi threads 10 (e -10); provide SRR ID dir that contains SRR_ID.sra
time fasterq-dump SRR11184871 -e 10 --split-files -O test_new_fastq_dir -o outfile_SRR11184871_new_test
#time fasterq-dump SRR17062757/SRR17062757.sra -e 10 --split-files -O SRR17062757_fastq_dir -o outfile_SRR17062757_test

#gzip fastq files
#gzip test_new_fastq_dir/*.fastq

#parallel implementation of gzip
pigz *fastq
pigz -p 8 ${sample_ID}-tmp/*.fastq


#confirm if the library is single or pair end:
esearch -db sra -query SRR11184871  | efetch -format runinfo | cut -d , -f 16

##confirm if the library is single or pair end: (both command works)
efetch -db sra -id SRR11184871 -format docsum | grep "LIBRARY_LAYOUT" -A 2 -m 2



########################
#SRA download
#You can use xargs and the sra-toolkit prefetch 
#to download every SRR id contained in a txt file list, like:
xargs -n1 prefetch < SRR_Acc_List.txt

#query ID
srx_id=SRX7805363

#query SRX ID to generate SRR ID for fetching SRA files
srr_id=$(esearch -db sra -query ${srx_id} | efetch -format runinfo| cut -d ',' -f 1 | grep SRR)

#prefetch saves time and fastens the process
time prefetch ${srr_id}
prefetch $(<SRA_Acc_List.txt) #List of IDs

#xargs -n1 prefetch < SRR_Acc_List.txt #from list

#check file integrity and check md5sums
vdb-validate ${srr_id}

#fasterq with multi threads 10 (e -10); provide SRR ID dir that contains SRR_ID.sra
time fasterq-dump ${srr_id} -e 10 --split-files -O test_new_fastq_dir -o outfile_SRR11184871_new_test

# #gzip fastq files
# gzip test_new_fastq_dir/*.fastq

#parallel implementation of gzip
pigz *fastq

#confirm if the library is single or pair end:
esearch -db sra -query ${srr_id}  | efetch -format runinfo | cut -d , -f 16

##confirm if the library is single or pair end: (both command works)
efetch -db sra -id ${srr_id} -format docsum | grep "LIBRARY_LAYOUT" -A 2 -m 2




# Throw your SRR numbers into a file called SRR_list.txt, one number per line.
# Then add this to a file called get_SRR_data.sh

#    #!/usr/bin/bash

#     fastq-dump --split-3 $1

# and run on the command line with:
cat SRR_list.txt | xargs -n 1 bash get_SRR_data.sh



#Make a list.txt file containing a single column of SRA numbers to download.
#then:

for i in $(cat list.txt); do 
	echo $i; date; 
	fasterq-dump -S $i; 
done

#It works well to use NCBI's web interface to find SRA samples of interest, download and open findings in Excel, then copy single column containing SRA numbers and paste into list.txt using document editor such as vim.

#After downloading including the "R" can be nice:

# for i in *_1.fastq; do 
# 	mv $i ${i%_1.fastq}_R1.fastq; 
# done

# for i in *_2.fastq; do 
# 	mv $i ${i%_2.fastq}_R2.fastq; 
# done

# and zip:
# pigz *fastq

#############################
# If you have a multi-core machine using pigz is much faster than traditional gzip.

# pigz, which stands for parallel implementation of gzip, is a fully functional replacement for gzip that exploits multiple processors and multiple cores to the hilt when compressing data. pigz was written by Mark Adler, and uses the zlib and pthread libraries.

# Pigz ca be used as a drop-in replacement for gzip. Note than only the compression can be parallelised, not the decompression.
################################
# Using pigz the command line becomes

mysqldump "$database_name" | pigz > $BACKUP_DIR/$database_name.sql.gz


Step 2: make the bash script
Create a bash with these lines:

#!/bin/bash

FILENAME="files2download.txt"

LINES=$(cat $FILENAME)

for LINE in $LINES; do
    echo $LINE
    parallel-fastq-dump -s $LINE -t 8 -O ./ --gzip 
done

#You can use parallel.

parallel -j 3 fastq-dump {} ::: SRR10611214 SRR10611215 SRR10611215 SRR10611216 SRR10611217
#The option -j says how many jobs should maximal run parallel. 
#So in this case maximal 3 identifier would be handled at the same time.

#How many jobs you can run parallel depends on your machine.

#You can also take a look at parallel-fastq-dump.


Converting to FASTQPermalink
I use the following bash script for converting multiple SRA to fastq. The script will search through my SRA directory for .SRA files and provide these files as the input for fasterq-dump

#Copy
#!/bin/bash
#fasterq-dump script
#loadmodules
module load sratoolkit/2.9.4
module load parallel
#find directory
find $SCRATCH/opt/sra/*sra | parallel 'fasterq-dump -O $SCRATCH/opt/fastq {}'



# (entrez-direct)[schhetri@ tools]$ time fasterq-dump SRR17062757/SRR17062757.sra -e 10 --split-files -O SRR17062757_fastq_dir -o outfile_SRR17062757_test
# spots read      : 25,243,584
# reads read      : 50,487,168
# reads written   : 50,478,374
# fasterq-dump SRR17062757/SRR17062757.sra -e 10 --split-files -O  -o   345.37s user 70.77s system 344% cpu 2:00.96 total
# (entrez-direct)[schhetri@ tools]$ ls
# SRR17062757            SRR17062757_fastq_dir  sratoolkit.3.0.1-mac64
# (entrez-direct)[schhetri@ tools]$ ls SRR17062757_fastq_dir
# outfile_SRR17062757_test         outfile_SRR17062757_test_1.fastq outfile_SRR17062757_test_2.fastq



Automating downloads using Python
#Since there are lots of SRA files associated with our samples, it would take a long time to manually run prefetch and fastq-dump for all the files. To automate this process, I wrote a small script in python to first download each SRA file using prefetch and then run fastq-dump. The code is shown below and also provided in this repo as fastq_download.py:

import subprocess

# samples correspond to Het_1, Het_2, Imm_1, Imm_2
sra_numbers = [
    "SRR2121685", "SRR2121686", "SRR2121687", "SRR2121688"
    ]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
