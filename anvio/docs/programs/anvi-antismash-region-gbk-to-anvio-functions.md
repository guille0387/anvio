This program **converts antiSMASH per-region GenBank files into a %(functions-txt)s for an anvi'o %(contigs-db)s.**

It is designed for the contigs-db analysis workflow where antiSMASH was run on GenBank records exported from anvi'o. The script reads only top-level files whose names end in `contig<number>.region<number>.gbk`.

The program consumes an %(antismash-region-gbk)s artifact and writes a standard five-column %(functions-txt)s that can be imported with %(anvi-import-functions)s.

## How it maps genes

For every CDS with an antiSMASH `/gene_kind`, the program reads `/locus_tag`. If the locus tag exactly matches `anvio_gene_<integer>`, only the integer is written to the `gene_callers_id` column.

For example:

{{ codestart }}
/locus_tag="anvio_gene_12401"
{{ codestop }}

becomes:

{{ codestart }}
12401
{{ codestop }}

This keeps the output compatible with the gene caller IDs already stored in the target %(contigs-db)s.

## Accession assignment

The `accession` column combines the region file stem and the antiSMASH `contig_edge` value.

For example:

{{ codestart }}
Day17a_QCcontig235.region001.gbk
{{ codestop }}

becomes the region id:

{{ codestart }}
Day17a_QCcontig235.region001
{{ codestop }}

and the final accession:

{{ codestart }}
Day17a_QCcontig235.region001_is_edge_TRUE
{{ codestop }}

If `contig_edge` is `False`, the accession ends in `_FALSE`. If the qualifier is missing or unrecognized, the script writes `_UNKNOWN`.

## Run the parser

Here is a typical run:

{{ codestart }}
python3 /Users/kpf734/antismash-anvio-parser/antismash_region_gbk_to_anvio_functions.py \
    -i /path/to/antismash-output-dir \
    -o antiSMASH.functions.tsv
{{ codestop }}

The output is a regular %(functions-txt)s with columns:

{{ codestart }}
gene_callers_id    source    accession    function    e_value
{{ codestop }}

The program writes two rows per importable CDS:

* `antiSMASH_type`, where the function is the antiSMASH region product.
* `antiSMASH_function`, where the function is the antiSMASH `gene_kind`.

## Import into anvi'o

Once the %(functions-txt)s is ready, import it into the matching %(contigs-db)s:

{{ codestart }}
anvi-import-functions \
    -c %(contigs-db)s \
    -i antiSMASH.functions.tsv
{{ codestop }}

The `e_value` column is blank because antiSMASH does not provide a real e-value for `/gene_kind` annotations.
