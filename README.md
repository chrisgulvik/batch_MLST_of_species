# batch MLST of species
__Purpose:__ to get allele (nucleotide) sequences according to multi-locus sequence types (STs)

## Input
- a tab-delimited file (despite its extension enigmatically being TXT) from [PubMLST](http://pubmlst.org/data/), labeled as a 'profile' for a particular species database
- all multi-FastA formatted files of alleles for the specified profile in the `cwd`; default file extension is TFA as [PubMLST](http://pubmlst.org/data/) uses, however any extension can be specified with the `--ext` opt

## Output
- multi-FastA format file of each allele for a ST  (ST-##.fas)
- a FastA file of the concatenated alleles for a ST  (concat_ST-##.fas)

#### Usage
    python get_allele_seqs_from_PubMLST_db -p profile.txt
