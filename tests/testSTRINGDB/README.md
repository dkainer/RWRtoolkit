Toy network from STRING-DB
==========================

1. Go to <string-db.org>.
2. Search for protein 'hexokinase' in organism 'homo sapiens'.
3. Settings > minimum required interaction score > highest confidence (0.900)
4. Click `+ More` button once (or twice ... don't recall for sure).
5. Exports > ... as simple tabular text output > download
6. Open the file `string_interactions.tsv` and separate the table into
   'layers' according by filtering each layer column (numerics). If the
   edge in that later has a score >0, keep the layer in that edge. Save
   these layers to separate files.
7. Make an `Rdata` file using `RWR_make_multiplex.R` and `flist.tsv`.
