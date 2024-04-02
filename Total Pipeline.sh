#!/bin/bash

#PBS -l select=1:ncpus=150:mem=500gb
#PBS -l walltime=20:00:00
#PBS -N rscript_job
module load anaconda3/personal
source activate myenv123





# Download samtools in the environment
conda install -n Yunkai bioconda::samtools -y
# Create a file called "BAM" to save the original 159 files
mkdir BAM
# Download the data from HPC to the local "BAM" folder
## Download data of patients 1 to 11


# Download QDNAseq file "QDNAseq_from_bam_chrX.R", shallowHRD file "312shallowHRD_hg19_1.13_QDNAseq_chrX.R", and "./cytoband_adapted_hg19.csv" file
BAM_BASE_DIR="./"
QDNA_SEQ_SCRIPT="./QDNAseq_from_bam_chrX.R"
SHALLOW_HRD_SCRIPT="./312shallowHRD_hg19_1.13_QDNAseq_chrX.R"
CYTOBAND_FILE="./cytoband_adapted_hg19.csv"


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
## Download "ks" with the version of 1.14.2
R -e "devtools::install_version('ks', version = '1.14.2')"





cd /rds/general/ephemeral/user/yt1823/ephemeral/pipeline1415
mkdir BAM


##copy files to BAM
cp /rds/general/project/pc_multiomics_project1/live/MRes/Yunkai/Data/all22/{14_T16-156_FT13.dups_marked.sorted.bam,14_T16-156_FT16.dups_marked.sorted.bam,14_T16-156_FT1.dups_marked.sorted.bam,14_T16-156_FT3.dups_marked.sorted.bam,14_T16-156_FT8.dups_marked.sorted.bam,14_T16-156_SN16-232_BC.dups_marked.sorted.bam,15_T16-186_BC1.dups_marked.sorted.bam,15_T16-186_FT1.dups_marked.sorted.bam,15_T16-186_FT2.dups_marked.sorted.bam,15_T16-186_FT5.dups_marked.sorted.bam,15_T16-186_FT6.dups_marked.sorted.bam,15_T16-186_FT8.dups_marked.sorted.bam} /rds/general/ephemeral/user/yt1823/ephemeral/pipeline1415/BAM/







# Download QDNAseq file "QDNAseq_from_bam_chrX.R", shallowHRD file "312shallowHRD_hg19_1.13_QDNAseq_chrX.R", and "./cytoband_adapted_hg19.csv" file
BAM_BASE_DIR="./"
QDNA_SEQ_SCRIPT="./QDNAseq_from_bam_chrX.R"
SHALLOW_HRD_SCRIPT="./312shallowHRD_hg19_1.13_QDNAseq_chrX.R"
CYTOBAND_FILE="./cytoband_adapted_hg19.csv"


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

# Downsampling 1-12 BAM files in the list
samples=$(sed -n '26,37p' bam_file_names.txt)



# Conduct parallel downsampling, the number of parallel processes can change, for exaple, 24 here
parallel -j 144 downsampling {1} {2} {3} ::: $samples ::: 0.6667 0.3333 0.1667 0.0333 0.0267 0.02 0.0167 0.0133 0.01 0.0067 0.0033 0.00167 :::+ 20x 10x 5x 1x 0.8x 0.6x 0.5x 0.4x 0.3x 0.2x 0.1x 0.05x





# QDNAseq 

# Defining the function for QDNAseq
processQDNAseq() {
    local bamFile=$1         
    local outputDir=$2       

    mkdir -p "$outputDir"    
    Rscript QDNAseq_from_bam_chrX.R "$bamFile" "$outputDir" 50
}

export -f processQDNAseq    # 导出函数以便parallel可以使用

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

# 生成任务列表
for downsampleDir in "${!downsampleToQDNA[@]}"; do
    outputDir="${downsampleToQDNA[$downsampleDir]}"  # 根据关联数组计算输出目录
    for bamFile in ./"$downsampleDir"/*sorted.bam; do
        # 写入任务到文件：BAM文件路径 和 输出目录
        echo "$bamFile $outputDir" >> "$taskFile"
    done
done

# Use parallel processing to run QDNAseq
parallel --colsep ' ' -j 144 processQDNAseq {1} {2} :::: "$taskFile"




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

