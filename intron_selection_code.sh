# Set the GTF input file used for alignment of human samples
gtf_input="Drosophila_melanogaster.BDGP6.90.gtf"


# Script to extract and sort genomic features (exons, introns, and UTRs) from a GTF file.
# Usage: ./bedtools.sh


# Define output paths for BED files by replacing the extension in the input filename
bed_gene="${gtf_input/.gtf/.gene.bed}"
bed_exon="${gtf_input/.gtf/.exon.bed}"
bed_5prime_utr="${gtf_input/.gtf/.5prime_utr.bed}"
bed_3prime_utr="${gtf_input/.gtf/.3prime_utr.bed}"
bed_intron="${gtf_input/.gtf/.intron.bed}"

# Function to process genomic features and output to BED files
process_features() {
    local feature_type=$1
    local feature_label=$2
    local output_file=$3

    echo "Processing $feature_type..."
    cat "$gtf_input" | awk -v label="$feature_label" 'BEGIN {OFS="\t"} $3 == label {print $1, $4, $5, label, $10, $7}' | bedtools sort > "$output_file"
}

# Extract and sort gene features (gene is labelled as intron  here; after proccessing only introns remain)
process_features "gene" "gene" "$bed_gene"

# Extract and sort exon features
process_features "exon" "exon" "$bed_exon"

# Extract and sort 5' UTR features
process_features "five_prime_utr" "5prime_utr" "$bed_5prime_utr"

# Extract and sort 3' UTR features
process_features "three_prime_utr" "3prime_utr" "$bed_3prime_utr"

# Function to subtract exons and UTRs from gene coordinates to get introns
subtract_features() {
    echo "Subtracting features to identify introns..."
    cat "$bed_gene" | bedtools subtract -s -a stdin -b "$bed_exon" | bedtools subtract -s -a stdin -b "$bed_5prime_utr" | bedtools subtract -s -a stdin -b "$bed_3prime_utr" | bedtools sort | bedtools merge -i stdin -s -c 4,5,6 -o distinct > "$bed_intron"
}

# Call function to process intron subtraction and merging
subtract_features

echo "Feature processing complete."
