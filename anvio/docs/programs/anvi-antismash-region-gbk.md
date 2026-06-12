This artifact is a collection of antiSMASH per-region GenBank files produced from contigs exported from anvi'o.

The parser reads only top-level files whose names end in `contig<number>.region<number>.gbk`, such as:

{{ codestart }}
Day17a_QCcontig235.region001.gbk
Day17a_QCcontig104.region001.gbk
Day17a_QCcontig74.region001.gbk
{{ codestop }}

The parser ignores combined antiSMASH GenBank files, JSON files, HTML files, JavaScript assets, and files in nested output directories.

For a CDS feature to become an anvi'o function annotation, it must include:

* a `/gene_kind` qualifier from antiSMASH, and
* a `/locus_tag` qualifier that starts with `anvio_gene_` and ends with an integer.

For example, this locus tag:

{{ codestart }}
/locus_tag="anvio_gene_12401"
{{ codestop }}

will be written to the output %(functions-txt)s as the anvi'o `gene_callers_id` value `12401`.
