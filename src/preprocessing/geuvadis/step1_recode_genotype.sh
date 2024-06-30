#!/usr/bin/env bash



# commonpath is the base directory path provided as a command-line argument
commonpath="$1"

# inputfilepath is the directory containing raw genotype VCF files
inputfilepath="${commonpath}raw/genotypes/"

# outputfilepath is the directory where the processed files will be saved
outputfilepath="${commonpath}preprocessed/recoded/"

# sampleinfofile contains sample information
sampleinfofile="${commonpath}raw/info/E-GEUV-1.sdrf.txt"

# idconvertfile is used to update SNP IDs in the VCF files
idconvertfile="${inputfilepath}dbSnp137_idconvert_unique.txt"

# BEST_SNP contains the list of SNPs to be retained in the analysis (best SNPs)
BEST_SNP="${inputfilepath}uniq_best_ids.txt"

# ALL_SNP contains the list of SNPs to be retained in the analysis (all SNPs)
ALL_SNP="${inputfilepath}uniq_all_ids.txt"

# header_file is the file containing header information for VCF files
header_file="${inputfilepath}header.txt"

if [ ! -d "$outputfilepath" ]; then
    mkdir "$outputfilepath"
fi

shopt -s nullglob

# Get the list of input files without sorting
input_files=("${inputfilepath}"*.vcf.gz)

for f in "${input_files[@]}"; do
    filename=$(basename "$f")
    outputname1="${filename%%.vcf*}_best"
    outputname2="${filename%%.vcf*}_all"

     PLINK 2 command to process VCF file
    # plink2: Command to run PLINK 2.0, a toolset for whole-genome association and population-based linkage analysis.
    # --vcf "$f": Specifies the input file format and file path. $f is the variable holding the path to the current VCF file to be processed.
    # --export Av: Specifies the output format. A means the data will be in allele format, and v indicates transposed data.
    # --update-name "$idconvertfile" 1 2: Updates the SNP IDs in the VCF file. The file specified by $idconvertfile contains two columns:
    #     Column 1: The current SNP ID.
    #     Column 2: The new SNP ID to replace the current one.
    # --extract "$BEST_SNP": Extracts a subset of SNPs. The file specified by $BEST_SNP contains the list of SNPs to be retained in the analysis.
    # --geno 0.05: Excludes SNPs with a missing genotype rate higher than 5%. SNPs missing in more than 5% of samples will be filtered out.
    # --out "${outputfilepath}${outputname1}": Specifies the base name for output files. The output files will be saved to the path and name given by combining $outputfilepath and $outputname1.

    plink2 --vcf "$f" --export Av --update-name "$idconvertfile" 1 2 --extract "$BEST_SNP" --geno 0.05 --out "${outputfilepath}${outputname1}"
    plink2 --vcf "$f" --export Av --update-name "$idconvertfile" 1 2 --extract "$ALL_SNP" --geno 0.05 --out "${outputfilepath}${outputname2}"
done
