# These pipeline was built for running in my environment, minor adjustment should be made when running in a different environment


# I run in PBS job submission system, I applied for the following amount of resources to process 10-20 WGS data (in BAM format) each time.
#!/bin/bash
#PBS -l select=1:ncpus=220:mem=800gb
#PBS -l walltime=20:00:00
#PBS -N rscript_job



# My working environment was myenv123
module load anaconda3/personal
source activate myenv123


# But if you want to use a new envrionment, for example, called "Yunkai", please download the following packages.

# Download samtools in the environment
conda install -n Yunkai bioconda::samtools -y

# Download the packages for QDNAseq
conda install -n Yunkai -y bioconda::bioconductor-biobase
conda install -n Yunkai -y bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg19
conda install -n Yunkai -y bioconda::bioconductor-qdnaseq
conda install -n Yunkai -y bioconda::bioconductor-dnacopy
conda install -n Yunkai -y bioconda::bioconductor-cghcall
conda install -n Yunkai -y bioconda::bioconductor-qdnaseq.hg19


# Download the packages for shallowHRD
conda install -n Yunkai -y conda-forge::r-ggpubr
conda install -n Yunkai -y bioconda::bioconductor-genomicranges
conda install -n Yunkai -y conda-forge::r-gridextra
conda install -n Yunkai -y conda-forge::r-desctools
conda install -n Yunkai -y conda-forge::r-ggrepel
## Download "devtools" package for downloading the "ks" package with the version of 1.14.2
conda install -n Yunkai -y conda-forge::r-devtools
# Download "ks" with the version of 1.14.2
R -e "devtools::install_version('ks', version = '1.14.2')"


# Put these files in your working directory: QDNAseq file "QDNAseq_from_bam_chrX.R", shallowHRD file "312shallowHRD_hg19_1.13_QDNAseq_chrX.R", and "./cytoband_adapted_hg19.csv" file 
BAM_BASE_DIR="./"
QDNA_SEQ_SCRIPT="./QDNAseq_from_bam_chrX.R"
SHALLOW_HRD_SCRIPT="./312shallowHRD_hg19_1.13_QDNAseq_chrX.R"
CYTOBAND_FILE="./cytoband_adapted_hg19.csv"
# Rememver also put the name list of the WGS data that you are working on. Here I put my list called "bam_file_names.txt"


# Create a file called "BAM" to save WGS BAM files. In this study, original 159 files were saved and used.
mkdir BAM

# Download the data from HPC to the local "BAM" folder. Here we put all 12 files from patient 6 and patient 7. You can adjust the paths as you want to test different WGS data.
cp /rds/general/project/pc_multiomics_project1/live/MRes/Yunkai/Data/all22/{6_T14-125_157_BC.dups_marked.sorted.bam,6_T14-125_62_FT6.dups_marked.sorted.bam,6_T14-125_63_FT7.dups_marked.sorted.bam,6_T14-125_64_FT8.dups_marked.sorted.bam,6_T14-125_65_FT9.dups_marked.sorted.bam,7_T14-137_158_BC.dups_marked.sorted.bam,7_T14-137_67_FT2.dups_marked.sorted.bam,7_T14-137_68_FT3.dups_marked.sorted.bam,7_T14-137_69_FT5.dups_marked.sorted.bam,7_T14-137_70_FT7.dups_marked.sorted.bam,7_T14-137_71_FT9.dups_marked.sorted.bam,7_T17-067_140_A.dups_marked.sorted.bam,7_T17-067_141_B.dups_marked.sorted.bam,7_T17-067_142_C1.dups_marked.sorted.bam,7_T17-067_143_D.dups_marked.sorted.bam,7_T17-067_144_E.dups_marked.sorted.bam,7_T17-067_145_F.dups_marked.sorted.bam,7_T17-067_146_G.dups_marked.sorted.bam} /rds/general/ephemeral/user/yt1823/ephemeral/pipeline67/BAM/




# Downsampling

downsampling(){
  sampleName=$1
  rate=$2
  key=$3  
  outputDir="./${key}BAM"
  mkdir -p "$outputDir"
  samplePath="./BAM/${sampleName}"
  outputFile="${outputDir}/$(basename "${sampleName}" .bam)_${rate}_${key}.bam"
  
  samtools view -b -s $rate "$samplePath" > "$outputFile"
  samtools sort "$outputFile" -o "${outputFile%.*}_sorted.bam"
  samtools index "${outputFile%.*}_sorted.bam"
  rm "$outputFile"
}
export -f downsampling

# Downsampling 124th-141st BAM files from the "bam_file_names.txt"
samples=$(sed -n '124,141p' bam_file_names.txt)



# Conduct parallel downsampling, the number of the parallel task running together each time can be changed, here I recommed to set the number of j to 11 x the number of WGS file processed, as there are 11 downsampling depths.
parallel -j 198 downsampling {1} {2} {3} ::: $samples ::: 0.6667 0.3333 0.1667 0.0333 0.0267 0.02 0.0167 0.0133 0.01 0.0067 0.0033 0.00167 :::+ 20x 10x 5x 1x 0.8x 0.6x 0.5x 0.4x 0.3x 0.2x 0.1x 0.05x





# QDNAseq 

# Defining the function for QDNAseq
processQDNAseq() {
    local bamFile=$1         
    local outputDir=$2       

    mkdir -p "$outputDir"    
    Rscript QDNAseq_from_bam_chrX.R "$bamFile" "$outputDir" 50
}

export -f processQDNAseq   
# Defining an associative array to map downsampled directories to QDNA output directories
declare -A downsampleToQDNA=(
    ["BAM"]="30xQDNA"
    ["20xBAM"]="20xQDNA"
    ["5xBAM"]="5xQDNA"
    ["1xBAM"]="1xQDNA"
    ["0.8xBAM"]="0.8xQDNA"
    ["0.6xBAM"]="0.6xQDNA"
    ["0.5xBAM"]="0.5xQDNA"
    ["0.4xBAM"]="0.4xQDNA"
    ["0.3xBAM"]="0.3xQDNA"
    ["0.2xBAM"]="0.2xQDNA"
    ["0.1xBAM"]="0.1xQDNA"
    ["0.05xBAM"]="0.05xQDNA"
)

# Prepare task file
taskFile="tasks.txt"
> "$taskFile"  

# Generate task file
for downsampleDir in "${!downsampleToQDNA[@]}"; do
    outputDir="${downsampleToQDNA[$downsampleDir]}" 
    for bamFile in ./"$downsampleDir"/*sorted.bam; do
        echo "$bamFile $outputDir" >> "$taskFile"
    done
done

# Use parallel processing to run QDNAseq
parallel --colsep ' ' -j 198 processQDNAseq {1} {2} :::: "$taskFile"




# shallowHRD
## I have tried this part on hpc, but the parallel running is a problem, I don't know why it can only run 2 QDNA files at same time

declare -A downsampleToQDNA=(
    ["BAM"]="30xQDNA"
    ["20xBAM"]="20xQDNA"
    ["5xBAM"]="5xQDNA"
    ["1xBAM"]="1xQDNA"
    ["0.8xBAM"]="0.8xQDNA"
    ["0.6xBAM"]="0.6xQDNA"
    ["0.5xBAM"]="0.5xQDNA"
    ["0.4xBAM"]="0.4xQDNA"
    ["0.3xBAM"]="0.3xQDNA"
    ["0.2xBAM"]="0.2xQDNA"
    ["0.1xBAM"]="0.1xQDNA"
    ["0.05xBAM"]="0.05xQDNA"
)

export SHALLOW_HRD_SCRIPT CYTOBAND_FILE

process_qdna() {
    qdnaFile=$1
    qdnaOutputDir=$2
    fileName=$(basename "$qdnaFile")
    subDirName="${fileName%.*}"
    shallowOutputDir="${qdnaOutputDir%QDNA}shallow"
    specificShallowOutputDir="$shallowOutputDir/$subDirName"
    mkdir -p "$specificShallowOutputDir"
    Rscript $SHALLOW_HRD_SCRIPT "$qdnaFile" "$specificShallowOutputDir" "$CYTOBAND_FILE"
}

export -f process_qdna

for downsampleDir in "${!downsampleToQDNA[@]}"; do
    qdnaOutputDir="${downsampleToQDNA[$downsampleDir]}"
    shallowOutputDir="${qdnaOutputDir%QDNA}shallow"
    mkdir -p "$shallowOutputDir"

    #both this two version run 2 files at same time, and I don't know why
    find "$qdnaOutputDir" -type f | parallel process_qdna {} "$qdnaOutputDir"
    ##find "$qdnaOutputDir" -type f | parallel -j 8 process_qdna {} "$qdnaOutputDir"

done
