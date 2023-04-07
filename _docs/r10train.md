# Training for R10 using Nanopolish

It was a pain to figure out the nanopolish training as it is not documented.
With Jared Simpon's help finally figured out.
Documenting incase it is useful to someone (very brifly though).

We sequenced the human methylated and nonmethylated dataset from https://zymoresearch.eu/pages/dna-methylation-standards. 
The nonmethylated data is used for training neuclerotide model. Both methylated and nonmethylated are needed for training cpg methylation.


We will be using BLOW5 files to train. Get the r10 branch from https://github.com/hiruna72/nanopolish/ (the commit used: ce29c2dadb245857ad4a86a33b60724d2f312ca4) and compile.

# Neucleotide model

1. The base-model provided by ONT for dna_r10.4.1_e8.2_400bps wa used as the base model (available: https://github.com/nanoporetech/kmer_models). The values there are normalised. We convert to nanopolish model format and assume the file name is r10.4.1_400bps.nucleotide.9mer.model:

```
#ont_model_name r10_180mv_450bps_9mer
#kit    r10_450bps
#strand template
#k      9
#alphabet       nucleotide
kmer	level_mean	level_stdv	sd_mean	sd_stdv	weight
AAAAAAAAA	54.2049064106214	3.0	1	1	1	1
AAAAAAAAC	58.590797271313	3.0	1	1	1	1
AAAAAAAAG	55.9520679998075	3.0	1	1	1	1
AAAAAAAAT	58.4335200897368	3.0	1	1	1	1
AAAAAAACA	63.6599592820753	3.0	1	1	1	1
```

Yes, we are doing dna_r10.4.1_e8.2_400bps, but in the model above it should be r10_180mv_450bps_9mer to trick nanopolish train to work. level_mean from ONT model are converted to pA by using `ont_level*23.027081308793+96.631070539403`. level_std is set to 3.0 for all kmers. Other columns such as sd_mena are unsused, so leave it at 1.


2. Assume we have the nonmethylated dataset in a [BLOW file](https://github.com/hasindu2008/slow5tools) called nonmeth_reads.blow5. Basecalled with super accuracy through [buttery-eel](https://github.com/Psy-Fer/buttery-eel) using Guppy:
```
buttery-eel  -g /install/ont-guppy-6.4.2/bin/  --config dna_r10.4.1_e8.2_400bps_sup.cfg  --device 'cuda:all' -i  nonmeth_reads.blow5 -q 10 -o  nonmeth.fastq --port 5555
```
Will create two files nonmeth.pass.fastq and nonmeth.fail.fastq. We will be using the passed reads.

3. Now let us use Minimap2 to map. Used the [hg38 with no alternate contigs](bit.ly/hg38noAlt) as reference:
```
minimap2 -ax map-ont /genome/hg38noAlt.fa --secondary=no -t20 nonmeth.pass.fastq | samtools sort - > nonmeth_pass.bam
samtools index  nonmeth_pass.bam
```
4. Before running nanopolish train, we need to replaces the few 'M' bases in the hg38noAlt.fa with 'N' bases. Otherwise Nanopolish will error out thinking it is a methylated symbol. Assume that patched reeferncee is hg38noAlt_M_replaced_N.fa. Now do the training.
```
slow5-nanopolish-train/nanopolish index  nometh.pass.fastq --slow5 nometh_reads.blow5
slow5-nanopolish-train/nanopolish train -r nometh.pass.fastq -g hg38noAlt_M_replaced_N.fa -b nometh_pass.bam  -t 72 --train-kmers=all --input-model-filename ./r10.4.1_400bps.nucleotide.9mer.model -d trained_nuc_models/
```

if training is successful, you should see files like r10_450bps.nucleotide.9mer.template.roundX.model under trained_nuc_models and should have something like:
```
#model_name	r10_450bps.nucleotide.9mer.template.model
#kit	r10_450bps
#strand	template
#alphabet	nucleotide
AAAAAAAAA	54.2003	1.99094	1	1
AAAAAAAAC	57.2213	2.59858	1	1
AAAAAAAAG	55.7136	2.52011	1	1
AAAAAAAAT	57.2297	2.71002	1	1
```

Also there will be a file like r10_450bps.nucleotide.9mer.template.round0.summmary.tsv. Make sure the "is_trained" coulumn is set to 1:
```
model_short_name	kmer	num_matches	num_skips	num_stays	num_events_for_training	was_trained	trained_level_mean	trained_level_stdv
r10_450bps.t	CAAAAAAAA	0	0	0	1000	1	56.42	1.82
r10_450bps.t	TTCAAAAAA	0	0	0	1000	1	60.69	2.46
r10_450bps.t	AATAAAAAA	0	0	0	1000	1	57.26	3.36
r10_450bps.t	CCTAAAAAA	0	0	0	1000	1	60.58	2.35
r10_450bps.t	AAAAAAAAA	0	0	0	1000	1	54.20	1.99
```
If they are 0, something is wrong. Debugging time - these issues may help https://github.com/jts/nanopolish/issues/825, https://github.com/jts/nanopolish/issues/761,  https://github.com/jts/nanopolish/issues/1059, https://github.com/jts/nanopolish/issues/1064.


5. Then I used the last model trained in the above as the input and trained 10 more rounds.
```
slow5-nanopolish-train/nanopolish train -r nometh.pass.fastq -g hg38noAlt_M_replaced_N.fa -b nometh_pass.bam -t 72 --train-kmers=all --input-model-filename trained_nuc_models/r10_450bps.nucleotide.9mer.template.round4.model -d 
_more/ --rounds 10
```

6. The created trained_nuc_models_more/r10_450bps.nucleotide.9mer.template.round9.model which looked like below is the final model I use:
```
#model_name	r10_450bps.nucleotide.9mer.template.model
#kit	r10_450bps
#strand	template
#alphabet	nucleotide
AAAAAAAAA	53.9989	3.83325	1	1
AAAAAAAAC	55.769	3.20809	1	1
AAAAAAAAG	56.4699	5.17555	1	1
AAAAAAAAT	56.385	3.47242	1	1
AAAAAAACA	59.1981	4.42718	1	1
AAAAAAACC	62.2545	5.09254	1	1
AAAAAAACG	60.252	4.37649	1	1
AAAAAAACT	62.0946	4.71852	1	1
AAAAAAAGA	89.6554	29.6328	1	1
AAAAAAAGC	63.4326	4.10237	1	1
AAAAAAAGG	59.7029	4.50327	1	1
AAAAAAAGT	62.8706	3.91294	1	1
```

# Methylation model

1. We need both nonmethylated (used for nucleotide training above) and methylated data for this. We already basecalled nonmethylated data above. Now let us basecall the methylated data.
```
buttery-eel  -g /install/ont-guppy-6.4.2/bin/  --config dna_r10.4.1_e8.2_400bps_sup.cfg  --device 'cuda:all' -i  meth_reads.blow5 -q 10 -o  meth.fastq --port 5555

```
2. Now combined passed reads of nonmethylated and methylated into one fastq.
```
cat nometh.pass.fastq meth.pass.fastq > positive_and_negative_pass.fastq

```
3. Merge the methylated and non methylated BLOW5 files:
```
slow5tools merge meth_reads.blow5 nometh_reads.blow5 -o merged.blow5
```

4. Now we need to create a methylated reference - one with all CGs changed into MG. For this we use the hg38noAlt_M_replaced_N.fa we created above.

```
# make a copy with chr renamed to chr_meth
sed 's/chr/chr_meth_/g' hg38noAlt_M_replaced_N.fa > hg38noAlt_chr_renamed_to_meth.fa
# use seqtk to create single lined FASTA and replace all CG sites with MG.  
seqtk  seq -l0 hg38noAlt_chr_renamed_to_meth.fa | sed 's/CG/MG/g'  > hg38noAlt_chr_renamed_to_meth_methylated_cpg.fa
```
M is the nanopolish symbolfor methylated. The boundaries (C'\n'G) will not be propoerly replaced if you use a multiline FASTA - that is why seqtk is used here. 

5. Map the methylated reads to the renamed reference hg38noAlt_chr_renamed_to_meth.fa above (not the hg38noAlt_chr_renamed_to_meth_methylated_cpg.fa)

```
minimap2 -ax map-ont hg38noAlt_chr_renamed_to_meth.fa--secondary=no -t20 meth.pass.fastq | samtools sort - > meth_pass.bam
samtools index meth_pass.bam
```

6. The the BAMs from nonmethylated (nometh_pass.bam we created before for nucleotide) and methylated (meth_pass.bam) together.
```
samtools merge merged.bam nometh_pass.bam meth_pass.bam 
samtools index merged.bam
```

7. Now create the hybrid reference to be later used with nanopolish.
```
# 
seqtk  seq -l0 hg38noAlt_M_replaced_N.fa > hg38noAlt_tmp.fa
cat hg38noAlt_tmp.fa hg38noAlt_chr_renamed_to_meth_methylated_cpg.fa > hg38noAlt_hybrid.fa
```

8. Now we need to create a base CPG model for nanopolish. What I did was, I took the r10_450bps.nucleotide.9mer.template.round9.model which was trained above, created a new model with 'M' in the alphabet and put the same levels and stds as for the corresponding 'C' counterparts. Assume this file name as r10.4.1_400bps.cpg.9mer.model.
```
#model_name r10_450bps.cpg.9mer.template.model
#kit    r10_450bps
#strand template
#k      9
#alphabet       cpg
kmer	level_mean	level_stdv	sd_mean	sd_stdv	weight
AAAAAAAAA	53.998900	3.833250	1	1	1	1
AAAAAAAAC	55.769000	3.208090	1	1	1	1
AAAAAAAAG	56.469900	5.175550	1	1	1	1
AAAAAAAAM	55.769000	3.208090	1	1	1	1
AAAAAAAAT	56.385000	3.472420	1	1	1	1
AAAAAAACA	59.198100	4.427180	1	1	1	1
AAAAAAACC	62.254500	5.092540	1	1	1	1
AAAAAAACG	60.252000	4.376490	1	1	1	1
AAAAAAACM	62.254500	5.092540	1	1	1	1
AAAAAAACT	62.094600	4.718520	1	1	1	1
AAAAAAAGA	89.655400	29.632800	1	1	1	1
AAAAAAAGC	63.432600	4.102370	1	1	1	1
AAAAAAAGG	59.702900	4.503270	1	1	1	1
AAAAAAAGM	63.432600	4.102370	1	1	1	1
AAAAAAAGT	62.870600	3.912940	1	1	1	1
AAAAAAAMA	59.198100	4.427180	1	1	1	1
AAAAAAAMC	62.254500	5.092540	1	1	1	1
AAAAAAAMG	60.252000	4.376490	1	1	1	1
AAAAAAAMM	62.254500	5.092540	1	1	1	1
AAAAAAAMT	62.094600	4.718520	1	1	1	1
AAAAAAATA	58.516100	4.245120	1	1	1	1
AAAAAAATC	60.291300	4.159900	1	1	1	1
AAAAAAATG	60.182800	4.553570	1	1	1	1
AAAAAAATM	60.291300	4.159900	1	1	1	1
AAAAAAATT	59.724500	3.398340	1	1	1	1
AAAAAACAA	101.404000	18.021500	1	1	1	1
AAAAAACAC	108.078000	12.647600	1	1	1	1
AAAAAACAG	100.971000	18.010800	1	1	1	1
AAAAAACAM	108.078000	12.647600	1	1	1	1
AAAAAACAT	107.755000	11.699300	1	1	1	1
AAAAAACCA	101.966000	20.410900	1	1	1	1
AAAAAACCC	109.189000	13.997700	1	1	1	1
```

9. Now Nanopolish methylation training time. 
```
slow5-nanopolish-train/nanopolish index  positive_and_negative_pass.fastq --slow5 merged.blow5
slow5-nanopolish-train/nanopolish train -r positive_and_negative_pass.fastq -g hg38noAlt_hybrid.fa -b merged.bam  -t 72 --train-kmers=all --input-model-filename r10.4.1_400bps.cpg.9mer.model -d trained_meth_models --rounds 10
```

10. The created  trained_meth_models/r10_450bps.cpg.9mer.template.round9.model is used as the final model, which looked like:
```
#model_name	r10_450bps.cpg.9mer.template.model
#kit	r10_450bps
#strand	template
#alphabet	cpg
AAAAAAAAA	53.7537	3.65572	1	1
AAAAAAAAC	56.2305	4.32346	1	1
AAAAAAAAG	57.3002	5.67074	1	1
AAAAAAAAM	62.9934	8.11058	1	1
AAAAAAAAT	56.7203	4.17047	1	1
AAAAAAACA	58.8781	4.15414	1	1
AAAAAAACC	62.5381	5.1778	1	1
AAAAAAACG	60.3428	5.31056	1	1
AAAAAAACM	63.1525	4.73173	1	1
AAAAAAACT	61.8409	4.19311	1	1
AAAAAAAGA	92.1803	31.218	1	1
AAAAAAAGC	63.5452	4.65351	1	1
AAAAAAAGG	61.4242	5.45205	1	1
AAAAAAAGM	68.0082	6.17807	1	1
AAAAAAAGT	62.7944	4.12666	1	1
AAAAAAAMA	59.1981	4.42718	1	1
AAAAAAAMC	62.2545	5.09254	1	1
AAAAAAAMG	61.1943	3.33436	1	1
AAAAAAAMM	62.2545	5.09254	1	1
AAAAAAAMT	62.0946	4.71852	1	1
AAAAAAATA	58.8587	4.67559	1	1
AAAAAAATC	60.3588	4.00821	1	1
AAAAAAATG	60.477	5.42577	1	1
AAAAAAATM	61.3608	3.74302	1	1
AAAAAAATT	59.6004	3.20899	1	1
AAAAAACAA	99.7551	17.5451	1	1
AAAAAACAC	107.298	12.0526	1	1
AAAAAACAG	101.281	17.5249	1	1
AAAAAACAM	109.249	13.1459	1	1
AAAAAACAT	107.25	11.8888	1	1
```

Make sure that r10_450bps.cpg.9mer.template.round9.summmary.tsv has the colum "was_trained" set to 1 for k-mers with CG and MG. If not debugging time. Note: I had to manually replace the inbuilt .inl model in nanopolish with values from the trained nucletide model (https://github.com/hiruna72/nanopolish/commit/5c7d01372abaef2aa2de9fde7227b1b09d34ea88). 
If you are using the latest r10 branch from https://github.com/hiruna72/nanopolish/ this is already done for you.
___

As we used minimal data to train this, the models arenot perfect. Some k-mers did not even have enough coverage to train well. But the good thing is, now anyone who has  better data can use the models that I trained as the base model and train futher.





