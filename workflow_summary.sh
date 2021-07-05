# Overall workflow
# Step 1: Base calling with guppy (see SLURM script, guppy_run_SLURM.sh)
singularity exec tools/guppy/guppy361.sif guppy_basecaller -i . -s ./out --kit SQK-RNA002 --flowcell FLO-FLG001 --calib_detect

# Step 2: Adapter trimming with PoreChop (see SLURM script, poreChop_run_SLURM.sh)

# Step 3: Alignment of RNA reads with minimap2 
for items in *.fastq; do minimap2 -ax splice -uf -k14 ref.fa  ${items} | samtools sort -o ${items%.fastq}.sam; done
for items in *.sam; do samtools view -S -bh ${items} > ${items%.sam}.bam; done
for items in *.bam; do samtools stats $items  > ${items%.bam}.stats.csv; done

# Step 4. Data analysis in R (All scripts were uploaded in file R)
# Step 5. For antisense hotspot analysis extract sense and antisense transcripts from bam files
# 5.1. Sort bam files with samtools sort
for items in *.bam; do samtools sort $items -o ${items%.bam}.sorted.bam; done
# 5.2. Index bam files with samtools index
for items in *BR?.sorted.bam; do samtools index $items; done

# see https://www.bioinformatics.recipes/recipe/view/antisense/#code
# 5.3. Extract all features from the postitive (+) and negative (-) strand
cat NN2_ENO_curated.gtf |  awk '$7 == "+" { print $0 }' > forward_features.gtf
cat NN2_ENO_curated.gtf |  awk '$7 == "-" { print $0 }' > reverse_features.gtf
cat SG17M_ENO_curated.gtf |  awk '$7 == "+" { print $0 }' > forward_features.gtf
cat SG17M_ENO_curated.gtf |  awk '$7 == "-" { print $0 }' > reverse_features.gtf

# 5.4.
for items in *.sorted.bam; do bedtools intersect -a $items -b forward_features.gtf > ${items%.bam}.forward_overlap.bam; done
for items in *.sorted.bam; do bedtools intersect -a $items -b reverse_features.gtf > ${items%.bam}.reverse_overlap.bam; done


# 5.5. extract antisense reads
for items in *.forward_overlap.bam; do samtools view -b -F 4 -f 16 $items > ${items%.forward_overlap.bam}.forward_antisense.bam; done
for items in *.reverse_overlap.bam; do samtools view -b -F 4 -f 16 $items > ${items%.reverse_overlap.bam}.reverse_antisense.bam; done
# 5.6. extract sense reads
for items in *.forward_overlap.bam; do samtools view -b -F 4 -F 16 $items > ${items%.forward_overlap.bam}.forward_sense.bam; done
for items in *.reverse_overlap.bam; do samtools view -b -F 4 -F 16 $items > ${items%.reverse_overlap.bam}.reverse_sense.bam; done


# 5.7. Convert from bam to bed
for items in *sense.bam; do bedtools bamtobed -i $items > ${items%.bam}.bed; done


# 6. Obtain genome coverage information
for items in *forward_antisense.bam; do bedtools genomecov -ibam $items -bg  > ${items%.forward_antisense.bam}.genomecov.forward_antisense.bed; done
for items in *reverse_antisense.bam; do bedtools genomecov -ibam $items -bg  > ${items%.reverse_antisense.bam}.genomecov.reverse_antisense.bed; done

