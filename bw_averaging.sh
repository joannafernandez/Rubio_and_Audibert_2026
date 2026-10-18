#!/bin/bash

# Exit on error
set -euo pipefail

source /home/bafj7/anaconda3/etc/profile.d/conda.sh
conda activate bioenv

# Set paths
wigToBigWig="/home/bafj7/Desktop/wigToBigWig"
chrom_sizes="/home/bafj7/Desktop/cuttag/hg38.chrom.sizes"
bigwig_dir='/media/bafj7/Jo3/MANUAL_LAMINB/bam/run1'  

# Go to the directory with .bigWig files
cd "$bigwig_dir" || exit

# Step 1: Find base names (remove _R#)
base_names=$(ls *.bw | sed -E 's/_R[0-9]+\.bw//' | sort -u)

# Step 2: Loop through each group
for base in $base_names; do
    echo "Processing $base..."

    files=(${base}_R*.bw)

    if [ "${#files[@]}" -lt 2 ]; then
        echo "Skipping $base: fewer than 2 files"
        continue
    fi

    wig_out="${base}_averaged_signal.wig"
    wig_filtered="${base}_filtered.wig"
    bw_out="${base}_averaged_signal.bigWig"

    # Average using wiggletools
    wiggletools mean "${files[@]}" > "$wig_out"

    # Remove chrEBV from the wig file
    grep -v '^chrEBV' "$wig_out" > "$wig_filtered"

    # Convert filtered wig to bigWig
    "$wigToBigWig" "$wig_filtered" "$chrom_sizes" "$bw_out"

    echo "✅ Saved: $bw_out"

    # Clean up intermediate wig files
    rm -f "$wig_out" "$wig_filtered"
done

echo "🎉 All done!"
