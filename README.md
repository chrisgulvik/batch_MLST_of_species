# batch MLST of species
__Purpose:__ to get allele (nucleotide) sequences according to multi-locus sequence types (STs)

## Input
- a tab-delimited file (despite its extension enigmatically being TXT) from [PubMLST](http://pubmlst.org/data/), labeled as a 'profile' for a particular species database
- all multi-FastA formatted files of alleles for the specified profile in the `cwd`; default file extension is TFA as [PubMLST](http://pubmlst.org/data/) uses, however any extension can be specified with the `--ext` opt

## Output
- multi-FastA format file of each allele for a ST  (ST-##.fas)
- a FastA file of the concatenated alleles for a ST  (concat_ST-##.fas)

## Usage
    python get_allele_seqs_from_PubMLST_db.py -p profile.txt

## Example Install
    cd $HOME
    pip install biopython numpy
    git clone https://github.com/chrisgulvik/batch_MLST_of_species.git
    echo 'export PATH="$PATH:$HOME/batch_MLST_of_species"' >> $HOME/.bash_profile
    chmod u+x ~/batch_MLST_of_species/get_allele_seqs_from_PubMLST_db.py


## Example Usage
Step 1 - Download profile and corresponding alleles:

    mkdir ~/Cdiff_MLST && cd ~/Cdiff_MLST
    for f in {profiles/cdifficile.txt,alleles/cdifficile/adk.tfa,alleles/cdifficile/atpA.tfa,alleles/cdifficile/dxr.tfa,alleles/cdifficile/glyA.tfa,alleles/cdifficile/recA.tfa,alleles/cdifficile/sodA.tfa,alleles/cdifficile/tpi.tfa};
        do wget http://pubmlst.org/data/"$f";
    done

Step 2 - Get all allele sequences for each ST

    get_allele_seqs_from_PubMLST_db.py -p cdifficile.txt
    ls concat_ST-*.fas | xargs cat > 351_Cdiff_concat_MLSTs.fa
    rm *.fas

Step 3 - Align and infer phylogeny

    muscle -in 351_Cdiff_concat_MLSTs.fa -out 351_Cdiff_concat_MLSTs.fa.aln
    python ~/scripts/fasta2phylip.py 351_Cdiff_concat_MLSTs.fa.aln 351_Cdiff_concat_MLSTs.aln.phy
    # collapse invariant sites
    raxmlHPC-PTHREADS -T 30 -f a -s 351_Cdiff_concat_MLSTs.aln.phy -n Cdiff_concat_MLSTs -p 654321 -x 54321 -N 100 -m GTRGAMMA
