# pipeline for antisense analyse need  bam files , database.csv 
# also require expressed_list.bed generate with pol II ChIP-seq

echo "usage: bash antisense_pipeline.sh expressed_list.bed treated_mark species"

# parameters
if [ ! -n "$3" ];then
    echo "need parameters"
    exit 0
fi
expressed=$1
treated=$2
species=$3

treat_bam=`awk -F "," -v vawk="$treated" '{if ($2==vawk) print $1".bam"}' database.csv`
awk '{print $4}' $expressed > expression_list.txt
# move dictionary to work place (code)
source ~/.bashrc
bam=`awk -F"," 'NR>1 {print $1".bam"}' database.csv`


# antisense annotation 
samtools merge -@ 16 treated.bam $treat_bam
samtools index -@ 16 treated.bam
bamCoverage -b treated.bam -o treated_fw.bw -bs 10 --normalizeUsing None --filterRNAstrand forward -p 16
bamCoverage -b treated.bam -o treated_rev.bw -bs 10 --normalizeUsing None --filterRNAstrand reverse -p 16
bigWigToBedGraph treated_fw.bw treated_fw.bg
bigWigToBedGraph treated_rev.bw treated_rev.bg


python ~/antisense_code/sliding_.py -f treated_fw.bg -r treated_rev.bg -sp $species >>python.report
python ~/antisense_code/merge_.py  -f $expressed -i result.bed -e expression_list.txt -c ~/reference/$species/sizes.genome >>python.report

rm temp*

# make merged 
cat ~/reference/$species/RefSeq_NMR.bed >>sas.bed
cat merge.bed >>sas.bed
python ~/python/bed/bed_to_gff_transcripts.py sas.bed AST_ST.gff
mv sas.bed AST_ST.bed

featureCounts -T 16 -p -O -s 2 -a AST_ST.gff -t transcript -g transcript_id -o featureCounts.txt $bam 2>>featurecount.report

conda activate R_4
Rscript ~/antisense_code/antisense.R $treated $species
conda activate base

cp Differential/list300.txt expression_list.txt Differential/quantile

python ~/python/chromatin_RNAseq/database_rep.py

for i in rep_*
do
    cd $i 
    bash ~/antisense_code/heatmap.sh ../Differential/quantile/ 
    bash ~/antisense_code/plotProfile.sh quantile
    cd ..
done