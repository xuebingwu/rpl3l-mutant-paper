# RPL3L patient RNA-seq analysis with kallisto
# input: a folder with fastq.gz files 
# output: gene count and tpm files for each sample

# kallisto index 
# MANE Select transcript sequences (v1.4) 
# supplemented with plasmid sequences used in generating AC16 cell lines
kallisto_index=/home/xw2629/genomes/indexes/kallisto/MANEv1.4_Plasmid

if [ ! -f $kallisto_index ]; then
    # download and preprocess MANE select transcript sequences using the following script:
    # /home/xw2629/genomes/sequences/MANE.py
    human_rna=/home/xw2629/genomes/sequences/MANE.GRCh38.v1.4.refseq_rna.fna.symbol.fa.gz
    # plasmid sequences used in generating AC16 cell lines
    plasmid_seq=/home/xw2629/genomes/sequences/plasmid-shRNA-RPL3L.fa
    # build index for kallisto
    cd /home/xw2629/genomes/indexes/kallisto
    kallisto index -i MANEv1.4_Plasmid $human_rna $plasmid_seq
fi

# input folder with fastq.gz files
# T189M/D308N patient: MM031, MM032 (technical replicate)
# Controls: 226dpb, 6mpb, 127dpb (downloaded from European Nucleotide Archive 
# (ENA project PRJEB26969, sample ID: 5818sTS, 5828sTS, and 5836sTS) )
fastqFolder=/media/rna/sequencing_data/RPL3L/rnaseq-patient-sample/fastq

# output folder
output=/media/rna/sequencing_data/RPL3L/rnaseq-patient-sample/kallisto-20250727-MANEv1.4_Plasmid

# sample names
samples="MM031 MM032 226dpb 6mpb 127dpb"
# fastq files for each sample should have the sample name in the file name
# e.g., *JY020*.fastq.gz
# for paired-end data:
# *JY020*_L00*_R1_*.fastq.gz for read1
# *JY020*_L00*_R2_*.fastq.gz for read2

# number of threads for kallisto
thread=40

mkdir $output
cd $output 

for sample in $samples 
do
    echo $sample
    
    # get fastq files for each sample: the order is R1 R2 if paired-end
    files=`ls $fastqFolder/*$sample*.fastq.gz`
    
    echo $files
    
    if [ sample == "MM031" ] || [ sample == "MM032" ]
    then
        # kallisto, paired-end data, strand specific
        # --bootstrap-samples=100 --pseudobam
        kallisto quant -t $thread -i $kallisto_index -o $sample --rf-stranded $files 
    else
        # single-end data for controls
        # --bootstrap-samples=100 --pseudobam
        kallisto quant -t $thread -i $kallisto_index -o $sample --rf-stranded --single -l 200 -s 20 $files 
    fi
done

# merge gene counts and tpms for all samples

# create header for gene count and tpm files
echo -n -e "Gene" > header
for sample in $samples
do
    echo -n -e "\t$sample"
done  >> header
echo "" >> header 

countfile=gene.count.tsv
tpmfile=gene.tpm.tsv

# removing old files
if [ -f $countfile ]; then
    rm $countfile
fi
if [ -f $tpmfile ]; then
    rm $tpmfile
fi

# initialize count and tpm files with the first column: gene names
# the first column of abundance.tsv except the first row
# the same order for all samples
for sample in $samples
do 
    echo $sample
    cut -f 1 $sample/abundance.tsv | tail -n +2 > $countfile 
    cp $countfile $tpmfile 
    break
done

# cut and paste count/tpm columns from each sample
# abundance.tsv has 5 columns: target_id, length, effective_length, est_counts, tpm
# we will extract the 4th and 5th columns for counts and tpm respectively
# and paste them to the countfile and tpmfile
# the first row of abundance.tsv is the header, we will skip it
for sample in $samples
do 
    echo $sample
    # count: column 4 of kallisto output
    cut -f 4 $output/$sample/abundance.tsv | tail -n +2 > a
    paste $countfile a > b
    mv b $countfile 
    # tpm: column 5
    cut -f 5 $output/$sample/abundance.tsv | tail -n +2 > a
    paste $tpmfile a > b
    mv b $tpmfile 
done

# add header to the count and tpm files
cp header header2
cat $countfile >> header
mv header $countfile
cat $tpmfile >> header2
mv header2 $tpmfile  
rm a

# the output file gene.count.tsv is used for DESeq2 analysis
#geneScore -i /home/xw2629/xuebing/RPL3L/code/deseq2-D308_vs_infant.tsv -skip 1 -cScore 3 -d /home/xw2629/genomes/geneSets/msigdb_v2024.1.Hs_GMTs/msigdb.v2024.1.Hs.symbols.gmt -o msigdb-deseq2-D308-vs-infant-ks-20250727 -p 1e-8 -plot
 


