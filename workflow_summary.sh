# Overall workflow
# Step 1: Base calling with guppy (see SLURM script, guppy_run_SLURM.sh)
singularity exec tools/guppy/guppy361.sif guppy_basecaller -i . -s ./out --kit SQK-RNA002 --flowcell FLO-MIN106 --calib_detect

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

# Extract sense and antisense transcripts, see
# see https://www.bioinformatics.recipes/recipe/view/antisense/#code
