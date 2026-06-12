This program **converts antiSMASH per-contig-region GenBank files into a %(functions-txt)s for an anvi'o %(contigs-db)s.**

It is designed for the contigs-db analysis workflow where antiSMASH was externally run on a GenBank file exported from %(anvi-export-genbank)s. The program reads only top-level files whose names end in `contig<number>.region<number>.gbk`.

It writes a standard five-column %(functions-txt)s that can be imported with %(anvi-import-functions)s.


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

If `contig_edge` is `False`, the accession ends in `_FALSE`. If the qualifier is missing or unrecognised, the program writes `_UNKNOWN`.

## Run the program

Here is a typical run:

{{ codestart }}
anvi-export-antismash \
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
