# Script to extract and sort genomic features (exons, introns, and UTRs) from a GTF file.

gtf_input=${1:-"/cellfile/datapublic/rmarel_1/Internship/Poll_II_spd/Data/Intron_selection/Drosophila_melanogaster.BDGP6.90.gtf"}
output_folder=${2:-$(dirname "$gtf_input")}
input_file_full_name=$(basename "$gtf_input")
input_file_name="${input_file_full_name%.*}"

head -n 10 "$gtf_input"

# Define output paths for BED files
bed_gene="${output_folder}/${input_file_name}.gene.bed"
bed_exon="${output_folder}/${input_file_name}.gene.bed"
bed_5prime_utr="${output_folder}/${input_file_name}.5prime_utr.bed"
bed_3prime_utr="${output_folder}/${input_file_name}.3prime_utr.bed"
bed_intron="${output_folder}/${input_file_name}.intron.bed"

# Function to process genomic features and output to BED files
process_features() {
    local feature_type=$1
    local feature_label=$2
    local output_file=$3

    echo "Processing $feature_type..."
    cat "$gtf_input" | awk -v label="$feature_label" 'BEGIN {OFS="\t"} $3 == label {print $1, $4-1, $5, label, $10, $7}' | bedtools sort > "$output_file"
}

# Extract and sort gene features
process_features "gene" "gene" "$bed_gene"

# Extract and sort exon features
process_features "exon" "exon" "$bed_exon"

# Extract and sort 5' UTR features
process_features "five_prime_utr" "five_prime_utr" "$bed_5prime_utr"

# Extract and sort 3' UTR features
process_features "three_prime_utr" "three_prime_utr" "$bed_3prime_utr"


# Function to subtract exons and UTRs from gene coordinates to get introns
subtract_features() {
    echo "Subtracting features to identify introns..."
    cat "$bed_gene" | bedtools subtract -s -a stdin -b "$bed_exon" | bedtools subtract -s -a stdin -b "$bed_5prime_utr" | bedtools subtract -s -a stdin -b "$bed_3prime_utr" | bedtools sort | bedtools merge -i stdin -s -c 4,5,6 -o distinct > introns_selected
}

# Call function to process intron subtraction and merging
awk -v OFS='\t' '{$4="intron"; print}' "$bed_gene" > temp && mv temp "$bed_intron"

echo "Feature processing complete."
