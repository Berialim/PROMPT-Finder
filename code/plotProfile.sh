echo "parameter 1 = bed_folder"

source ~/.bashrc
# step 1 generate bw file and normalize them
mkdir bw
# generate none normalized bw file
declare -A reads
reads=()
for i in *.bam
do
    samtools index -@ 16 $i
    bamCoverage  -b $i  -o bw/${i%.*}_fw.bw -bs 10 --normalizeUsing None --filterRNAstrand forward -p 16 --ignoreDuplicates --samFlagExclude 260
    bamCoverage  -b $i -o bw/${i%.*}_rev.bw -bs 10 --normalizeUsing None --filterRNAstrand reverse -p 16 --ignoreDuplicates --samFlagExclude 260
    # calculate the total reads wiht samtools
    total_reads=`samtools view -@ 16 -F 260 -c $i`
    reads[$i]=$total_reads
done
cd bw
 

# normalize bw file with RPKM
for i in *_fw.bw
do
    python ~/python/deeptools/normalize_bw.py $i ${i%_*}_rev.bw ${reads[${i%_*}.bam]} 10 no_neg
done

echo "bw file done"
cd ..

# step 2 generate strand separate bed file
rm -r temp_fold
cp -r $1 temp_fold
cd temp_fold
for i in *.bed
do
    python ~/python/bed/remove_chrY.py $i
    python ~/python/bed/select_bed_strand.py $i
done
cd ..


# step 3 compute matrix and plot
mkdir plot
name1=`ls temp_fold | grep "_fw.bed" | sed "s/_fw.bed//g"`

for i in $name1
do
    bed="temp_fold/$i"
    name=$i
    fw=`ls bw | grep "fw.bw" | sed "s/.*_fw.bw/bw\/&/"`
    rev=`ls bw | grep "rev.bw" | sed "s/.*_rev.bw/bw\/&/"`
    fw_sample=`echo $fw | sed "s/\.bw//g" | sed "s/bw\///g"` 
    rev_sample=`echo $rev | sed "s/\.bw//g"  | sed "s/bw\///g"`
    # calculate sense counts in geneBody
    computeMatrix scale-regions -p 16 -b 5000 --regionBodyLength 10000  -S $fw --missingDataAsZero \
    -R ${bed}_fw.bed  \
    -o result_gene.fw.gz -p 16 
    computeMatrix scale-regions -p 16 -b 5000 --regionBodyLength 10000  -S $rev --missingDataAsZero \
    -R ${bed}_rev.bed  \
    -o result_gene.rev.gz -p 16 
    # merge sense 
    computeMatrixOperations rbind -m result_gene.fw.gz result_gene.rev.gz -o merged.gz 
    # change the labels
    computeMatrixOperations relabel -m merged.gz -o merged_gene.gz --sampleLabels $fw_sample
    rm merged.gz
    # calculate antisense counts upstream TSS
    computeMatrix scale-regions -p 16 -S $rev \
    -R ${bed}_fw.bed -b 5000 --regionBodyLength 10000 --missingDataAsZero \
    -o result.rev.gz -p 16 
    computeMatrix scale-regions -p 16 -S $fw \
    -R ${bed}_rev.bed -b 5000 --regionBodyLength 10000 --missingDataAsZero \
    -o result.fw.gz -p 16 
    # merge antisense gz
    computeMatrixOperations rbind -m result.rev.gz result.fw.gz -o merged.gz 
    # change labels
    computeMatrixOperations relabel -m merged.gz -o merged_antisense.gz --sampleLabels $rev_sample
    # reverse antisense gz file
    python ~/python/deeptools/Matrix_as_negative.py merged_antisense.gz
    rm merged_antisense
    computeMatrixOperations cbind -m merged_gene.gz merged_antisense.gz -o plot/${name}.gz
    # python ~/python/deeptools/remove_na_matrix.py plot/${name}.gz
    rm *.gz
done

cd plot
for i in *.gz
do
    computeMatrixOperations filterValues -m $i -o ${i%.*}_filterd.gz --max 5000 --min -5000 &
done
wait
for i in *.gz
do
    plotProfile -m $i -out ${i%.*}.pdf --plotFileFormat pdf --perGroup&
done
wait

echo "all finished"
