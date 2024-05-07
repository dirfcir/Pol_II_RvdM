from pathlib import Path

import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm
from collections import Counter

input_file_path = Path('Drosophila_melanogaster.BDGP6.90.gtf')

# %%
# We can parse the .gtf file as a normal tab-delim file for a quick inspection
df = pd.read_csv(input_file_path,
                 header=4,
                 names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                 delimiter="\t")
print(f"Num. genes: {sum(df['feature'] == 'gene')}")
print(f"Possible features: {set(df['feature'])}")

# %%
# Extracts genes, exons, 3'UTR and 5'UTR from .gtf file
gtf_records = BedTool(input_file_path)

genes = BedTool([record for record in tqdm(gtf_records, desc='Extracting genes') if record.fields[2] == 'gene']).sort()
exons = BedTool([record for record in tqdm(gtf_records, desc='Extracting exons') if record.fields[2] == 'exon']).sort()
three_prime_utr = BedTool(
    [record for record in tqdm(gtf_records, desc="Extracting 3' UTR") if record.fields[2] == 'three_prime_utr']).sort()
five_prime_utr = BedTool(
    [record for record in tqdm(gtf_records, desc="Extracting 5' UTR") if record.fields[2] == 'five_prime_utr']).sort()

# %%
# Compare to bash script results:
bash_script_genes = BedTool('Drosophila_melanogaster.BDGP6.90.gene.bed')
bash_script_exons = BedTool('Drosophila_melanogaster.BDGP6.90.exon.bed')
bash_script_3_prime_utr = BedTool('Drosophila_melanogaster.BDGP6.90.3prime_utr.bed')
bash_script_5_prime_utr = BedTool('Drosophila_melanogaster.BDGP6.90.5prime_utr.bed')

# This would fail, since the bash script is not keeping all the attributes (e.g. gene_name)
# assert bash_script_genes == genes

assert len(bash_script_genes) == len(genes)

merged_exons = exons.merge()
merged_bash_script_exons = bash_script_exons.merge()

for exon_1, exon_2 in zip(exons, bash_script_exons):
    assert exon_1.start == exon_2.start, f"{exon_1.start=} while {exon_2.start=}"
# %%

introns = genes.subtract(exons, s=True).subtract(three_prime_utr, s=True).subtract(five_prime_utr,
                                                                                   s=True).sort().merge()
introns.saveas('introns.bed')

introns_bash_script = BedTool("Drosophila_melanogaster.BDGP6.90.intron.bed")
# %%
lengths_bash = [len(x) for x in introns_bash_script]
lengths = [len(x) for x in introns]

counter_bash = dict(Counter(lengths_bash))
counter = dict(Counter(lengths))
# %%
# with open('/home/jakub/Desktop/Pol_II_RvdM/Antonis_code/Drosophila_melanogaster.BDGP6.90.intron_junction.bed') as antonis_file:
#     antonis_text = antonis_file.read()
    
# with open('/home/jakub/Desktop/Pol_II_RvdM/Drosophila_melanogaster.BDGP6.90.intron.bed') as rickfrid_file:
#     rickfrid_text = rickfrid_file.read()
# # %%

# df_antonis = pd.read_csv('/home/jakub/Desktop/Pol_II_RvdM/Antonis_code/Drosophila_melanogaster.BDGP6.90.intron_junction.bed', delimiter='\t', names=[str(x) for x in range(6)])

# df_ricfrid = pd.read_csv('/home/jakub/Desktop/Pol_II_RvdM/Drosophila_melanogaster.BDGP6.90.intron.bed', delimiter='\t', names=[str(x) for x in range(6)])

# df_antonis.iloc[0] == df_ricfrid.iloc[0]

# for index, row_antonis in df_antonis.iterrows():
#     row_ricfrid = df_ricfrid.loc[index]
#     if not all(row_ricfrid==row_antonis):
#         print("Mismatch")
#         break
    
    
    


    