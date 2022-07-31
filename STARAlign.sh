#Using STAR (Spliced Transcripts Alignment to a Reference) Tool to map the reads from a sequence library to a reference genome
#First create your Genome Index using the reference genome and the Gene Annotation File (GTF/GFF3)
# Note that Genome Indexes only to be generated once for every genome (e.g., Human, Mouse, Fruit Fly, etc.)
#Download reference genome and annotation file from the Ensembl website 
#Ensembl FTP Website: https://useast.ensembl.org/info/data/ftp/index.html
#To download a file: wget *link address*

#Download Miniconda if you dont have it already, if you do activate conda
source ~/.bashrc
#(base) now precedes the -bash-4.2 prompt
#If source ~/.bashrc does not run Conda, use: export PATH="/network/rit/lab/bioinformaticslab/2021_Bioinformatics_Program/YourFolderInTheBioinformaticsLabFolder/miniconda3/bin:$PATH"

#Download STAR
conda install -c bioconda star
#It will ask to download the following programs, type either yes or y and hit enter

#make a new directory (mkdir) in your folder to contain the output files from STAR
mkdir STAR

#Create STAR Alignment script
touch STARAlign.txt
nano STARAlign.txt

#specify the script (Below is what the script contained- For further explaination of the options used in the alignment script, refer to personal notes or STAR manual)
#Link to STAR manual (https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
#!/bin/bash
#SBATCH --job-name=STAR_Align      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --mem=50gb                     # Job memory request
#SBATCH --cpus-per-task=5            # Number of CPU cores per task
#SBATCH --time=36:00:00              # Time limit hrs:min:sec
#SBATCH --output=STARALIGN.log     # Standard output and error log

#source ~/.bashrc
#for i in $(ls /network/rit/lab/bioinformaticslab/Project1/ | grep '.fastq.gz' | sed s/[12].fastq.gz// | sort -u)
#do
 STAR --genomeDir /network/rit/lab/bioinformaticslab/HumanGenome  --runThreadN 5 --readFilesIn /network/rit/lab/bioinformaticslab/Project1/"$i"1.fastq.gz /network/rit/lab/bioinformaticslab/Project1/"$i"2.fastq.gz --outFileNamePrefix "$i" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand gunzip -c

#done

#STARAlign.txt (END)

#Next, run the alignment 
sbatch STARAlign.txt
#You can monitor your job on the server using squeue -u *Your NETID*

#After job has run, you can use the STARALIGN.log file to check to make sure it finished running successfully for all samples

