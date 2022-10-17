For Tidalbase TEs, non-*Dmel* seqs can be removed with the `seqkit` command below. 
This works for Tidalbase because *Dmel* TEs don't have 'Dmel' in the name, but other TEs have the
4 letter species abbreviation.

```bash
# seqkit 2.3.1
seqkit grep -r -n -v -p "\\|D.{3}\\\\" Tidalbase_transposon_sequence.fasta > repeats.fasta
```

After running the nextflow script, `fix_seqnames.R` can be run and used for `ptera`