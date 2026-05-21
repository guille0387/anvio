This program **converts antiSMASH output files (per-contig-region GenBank files) into a %(functions-txt)s**

It is designed for the contigs-db analysis workflow where antiSMASH was run externally on GenBank file exported from %(anvi-export-genebank), and extract key information from antiSMASH contig-region files.

The program consumes an %(antismash-contig-region-gbk)s artifact and writes a standard five-column %(functions-txt)s that can be imported with %(anvi-import-functions)s.

## Run the parser

Here is a typical run:

{{ codestart }}
python3 antismash_region_gbk_to_anvio_functions.py \
       -i /path/to/antismash-output-dir \
       -o antiSMASH.functions.tsv
{{ codestop }}

The output is a regular %(functions-txt)s with columns:

{{ codestart }}
gene_callers_id    source    accession    function    e_value
{{ codestop }}

The program writes two rows per importable CDS:

* `antiSMASH_type`, where the function is the antiSMASH region type.
* `antiSMASH_function`, where the function is the antiSMASH function type `gene_kind`.

The accession column is assigned from the antiSMASH contig-region id and the region’s contig_edge is TRUE or FALSE.


## Import into anvi'o

Once the %(functions-txt)s is ready, import it into the matching %(contigs-db)s:

{{ codestart }}
anvi-import-functions \
    -c %(contigs-db)s \
    -i antiSMASH.functions.txt other_functions.txt ... ...
{{ codestop }}

The `e_value` column is blank because antiSMASH does not provide a real e-value for `/gene_kind` annotations.
