#!/bin/bash

#Step 1 in this process is moving your downloaded FASTQ files into a single folder without subdirectories. This can be done manually by simply going to the download folder, searching for "fastq.gz" files in the search bar (ensuring that it is searching within your folder), and then copying those into another folder within your working directory. Alternatively, within the terminal window, you can query and move all FASTQs with a wildcard expression such as "cp */*.fastq.gz 'your directory'". Meaning, move all FASTQ.gz files within all subdirectories into the desination folder specified. 
#BEFORE DOING ANYTHING ELSE, cd into the folder you will be working in.
#This command activates the Qiime2 environment. Depending on the computer and installation, the environment name (in this case, “qiime2-2019.7”) may differ. 
source activate qiime2-2021.4;

#Here this is setting the taxonomy database variables as the path to their location on the computer. This will change from computer to computer but will be static once this path is established. 
gg='Path/To/GreenGenes/Database';
sv='Path/To/Silva/Database';

#These steps will prompt the user to input their files or variables after the script is launched. 
read -p "Input directory path where FASTQ files are contained (no subdirectories): " directory;
read -p "Input file path for map file (.txt): " mapfile;
read -p "Specify sampling depth for alpha and beta diversity statistics: " depth;

#Here this will reiterate everything you've entered and then begin the script. This is a good time to verify everything is correct before leaving the computer to run. 
echo "FASTQ directory is $directory";
echo "Map file path is $mapfile";
echo "Sampling depth is $depth";
echo "Database is Silva and Greengenes";

#This command is the start of the official workflow. We are importing our raw FASTQ read1 and read2 files downloaded from basespace here. The –-type and –-input-format denoted is specific to Illumina files and to paired end sequencing (this would be different for single end read data). The –-input-path you feed will be pulled from what was prompted above but it's just where all your FASTQ files can be found. The output is a .qza file which is specific to Qiime2 and cannot be directly used in other workflows like a .fna file might be (though export and conversion to biom/tsv files is simple).

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $directory  --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza;
echo "Sequences imported";

#This command is summarizing your demuxed read files. It will output reads for each sample and show you a quality map for all the reads. This can be useful for quickly checking how your samples did on the sequencer and can inform your choices on parameters in the next command (for example, truncating at a certain position due to lower quality scores). 
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-paired-end.qzv;

#The next command uses a dada2 plugin to filter out poor quality reads, determine the error rate and correct errors from the specific Illumina run, and remove chimeras. Additionally, it will dereplicate sequences which helps reduce file size in the long run. You should note that if you are processing data from multiple sequencing runs, each run should be run through this step separately and then combined together after. This is because each run will have its own specific error rate and when combined, data can be lost or kept when it should not have been. The input is your file from above. The four parameters directly following the input are to inform dada2 on where to trim or truncate your sequences so they are all the same size. The “trim” parameter is where you want to cut on the left side and “trunc” is where you want to cut on the right. For this run, we would like to keep the full-length read so we pick “0” and “250” as cutting locations. If you had primer sequences included in the read, you would need to use the "trim" parameter to remove these from each read (this is the case for Shoreline Biome reads for example).  

#IMPORTANT: You will want to assign the number of processors that will be working on this job. The “-p--threads” command will let you choose how many processors. Typing “0” will, counter intuitively, assign all available processors to the job. I would recommend you do this to maximize the speed. The “–-verbose” flag will narrate what is happening with the process and should let you monitor things so that you will be able to judge whether the process is hung up on something or not.

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 250 --p-trunc-len-r 250 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads 0 --verbose;
echo "dada2 filter and denoise complete";

#The next three commands generate visualization for your different objects. The first will show you your samples following filtering and denoising (how many reads remain). The second will show you your different replicated sequences (less useful), and the third will show you the denoising process including the number of reads lost during filtering, denoising, read merging, and chimera removal. 
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file $mapfile;
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv;
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization stats-dada2.qzv;

echo "Table and Rep Seqs Generated";

#This command will align all of your sequences using MAFFT, mask uninformative positions, generate an unrooted tree, and then finally root it. The output from this step is used for phylogentic alpha and beta diversity statistics. 
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --p-n-threads 10 --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#This command will generate a large variety of alpha and beta diversity statistics, as well as a rarefied table of your samples. Make sure that you enter the correct sampling depth based on the previous step. You will have entered this number in the script at the beginning and should have been informed by what you know about your sample type. You can always regenerate this with a different sampling depth later. All of your outputs will be put into the created directory. If you want to generate only specific metrics or a specific kind of metric, look on the Qiime2 tutorial pages for the desired command. 
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth $depth --m-metadata-file $mapfile --output-dir core-metrics-results;

#This next section will detail how you assign and visualize taxonomy of your data. This is done using your rep-seqs file generated from the dada2 step and a classifier. This classifier is basically a big library of sequences with a taxonomic ID assigned. This is a file you can download from the Qiime2 website or generate yourself. This script will generate Silva and GreenGenes taxonomy files and barplots by default. The output from these steps is a taxonomy file which shows the taxonomic assignment for each ASV. This tends to take the longest of all the commands. 
echo "Classification step";
qiime feature-classifier classify-sklearn --i-classifier $gg --i-reads rep-seqs.qza --p-n-jobs -1 --o-classification taxonomy-gg.qza;
qiime metadata tabulate --m-input-file taxonomy-gg.qza --o-visualization taxonomy-gg.qzv;
qiime feature-classifier classify-sklearn --i-classifier $sv --i-reads rep-seqs.qza --p-n-jobs -1 --o-classification taxonomy-sv.qza;
qiime metadata tabulate --m-input-file taxonomy-sv.qza --o-visualization taxonomy-sv.qzv;

#Here we take our generated taxonomy data and combine that with our original table and our metadata file to generate an interactive taxa bar plot. 

qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy-gg.qza --m-metadata-file $mapfile --o-visualization taxa-bar-plots-gg.qzv;
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy-sv.qza --m-metadata-file $mapfile --o-visualization taxa-bar-plots-sv.qzv;
#The interactive bar plots can be downloaded as SVG files or a CSV file can be generated with the raw numbers at a specific taxonomic level. 

#Import Qiime2 artifacts into R using the "Qiime2 to decontam, Phyloseq, DESeq2, pheatmap" walkthrough document. This imports your Qiime2 artifacts into phyloseq which can then be used in other R packages.

echo "Script complete";
#You can then go and view the plots generated from this. This is done by either typing "qiime tools view 'your .qzv file'" or going to https://view.qiime2.org/ and dragging your object into the window. Make sure you specify .qzv files rather than .qza files. The desired ones will also be labeled as “emperor” plots for the core metrics output. 

