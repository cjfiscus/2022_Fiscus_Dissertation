#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --output=./std/%j.stdout
#SBATCH --error=./std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="mk_data"
#SBATCH -p batch


# bedtools 2.28.0
module load bedops/2.4.24

ARAPORT11=/rhome/cfisc004/bigdata/data/Araport11
GENOME=/rhome/cfisc004/bigdata/data/genomes/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
SEQDIR=/rhome/cfisc004/bigdata/projects/kmers_arabidopsis/data/sequences

# prepare geneset 

## transposable element genes 
sed 's/Chr//g' "$ARAPORT11"/Araport11_transposable_element_gene.201606.bed > "$SEQDIR"/Araport11_transposable_element_gene.201606.bed
bedtools getfasta -name -fi "$GENOME" -bed "$SEQDIR"/Araport11_transposable_element_gene.201606.bed > "$SEQDIR"/Araport11_transposable_element_gene.201606.fa
rm "$SEQDIR"/Araport11_transposable_element_gene.201606.bed

## pseudogenes
grep "pseudogenic_transcript" "$ARAPORT11"/Araport11_GFF3_genes_transposons.201606.gff | sed 's/Chr//g' > "$SEQDIR"/Araport11_pseudogenes.gff
gff2bed < "$SEQDIR"/Araport11_pseudogenes.gff > "$SEQDIR"/Araport11_pseudogenes.bed
bedtools getfasta -name -fi "$GENOME" -bed "$SEQDIR"/Araport11_pseudogenes.bed > "$SEQDIR"/Araport11_pseudogenes.fa
rm "$SEQDIR"/Araport11_pseudogenes.gff "$SEQDIR"/Araport11_pseudogenes.bed

# cat geneset
cat "$ARAPORT11"/Araport11_genes.201606.cds.repr.fasta "$ARAPORT11"/Araport11_non_coding.201606.cdna.fasta \
	"$SEQDIR"/Araport11_transposable_element_gene.201606.fa "$SEQDIR"/Araport11_pseudogenes.fa \
	> "$SEQDIR"/Araport11_genes_all.fa 

# make gene list
rm "$SEQDIR"/gene_lst.txt # remove if exists
echo -e "id\tclass" >> "$SEQDIR"/gene_lst.txt
grep ">" "$ARAPORT11"/Araport11_genes.201606.cds.repr.fasta | cut -d"|" -f1 | sed 's/>//g' | sed 's/$/\tprotein-coding/g' >> "$SEQDIR"/gene_lst.txt
grep ">" "$ARAPORT11"/Araport11_non_coding.201606.cdna.fasta | sed 's/>//g' | sed 's/$/\tnon-coding/g' >> "$SEQDIR"/gene_lst.txt
grep ">" "$SEQDIR"/Araport11_transposable_element_gene.201606.fa | sed 's/>//g' | sed 's/$/\ttransposable element/g' >> "$SEQDIR"/gene_lst.txt
grep ">" "$SEQDIR"/Araport11_pseudogenes.fa | sed 's/>//g' | sed 's/$/\tpseudogene/g' >> "$SEQDIR"/gene_lst.txt

# prepare repeat library
cat "$SEQDIR"/athrep.ref "$SEQDIR"/5SrRNA.fasta "$SEQDIR"/45S_At3.fasta "$SEQDIR"/CEN_cluster_seqs.fa \
	> "$SEQDIR"/A_thal_repeat_lib.fa 

# cat seqs together
cat "$SEQDIR"/Araport11_genes_all.fa "$SEQDIR"/A_thal_repeat_lib.fa > "$SEQDIR"/A_thal_seqs_all.fa 
