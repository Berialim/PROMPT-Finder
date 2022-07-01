echo $1
for i in *.bam
do
    samtools index -@ 16 $i
    bamCoverage -b ${i} -o ${i%.*}.bw --binSize 40 \
    --normalizeUsing RPKM -p 16 --smoothLength 120
done
# plot heatmap for all coding genes and 
mkdir visualization
bw=`ls *.bw`

# comput matrix with quantile
cp -r $1 ./quantile
cd quantile
for i in *.txt
do
    python ~/python/bed/filter_bed.py $i ~/reference/hg38/RefSeq.bed inin ${i%.*}.bed
    python ~/python/bed/remove_chrY.py ${i%.*}.bed
done

cd ..
for i in quantile/*.bed
do
    ii=`echo $i|sed "s/.*\///g"`
    ii=`echo ${ii%.*}`
    # computeMatrix scale-regions -S $bw -R $i -b 5000 --regionBodyLength 10000 -a 5000 -o visualization/${ii}_ab.gz \
    # -p 16 2>>matrix.report 
    computeMatrix scale-regions -S $bw -R $i -b 5000 --regionBodyLength 10000 -o visualization/${ii}.gz \
    -p 16 2>>matrix.report 
done


cd visualization
for i in *.gz 
do
    plotHeatmap -m $i -out ${i%.*}_heatmap.pdf --colorMap Reds --whatToShow "heatmap and colorbar" --plotFileFormat pdf &
done
wait
for i in *.gz
do
    computeMatrixOperations filterValues -m $i -o ${i%.*}_filterd.gz --max 5000 &
done
wait
for i in *.gz
do
    plotProfile -m $i -out ${i%.*}_profile.pdf --perGroup  --plotFileFormat pdf &
    # plotProfile -m $i -out ${i%.*}.pdf --numPlotsPerRow 1  --plotFileFormat pdf &
done
wait

echo "all done"