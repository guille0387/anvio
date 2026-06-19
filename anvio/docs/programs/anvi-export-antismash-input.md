This program exports an annotated %(genbank-file)s from a %(contigs-db)s. It is intended for antiSMASH input, and it writes CDS features for coding gene calls in the contigs database.

To run it, provide an input database and an output path:

{{ codestart }}
anvi-export-antismash-input -c %(contigs-db)s \
                            -o path/to/output.gbk
{{ codestop }}

By default, all coding gene calls with sequences available in the contigs database are exported. You can limit the output to a subset of gene callers with `--gene-caller`, and you can skip shorter contigs with `--min-contig-length`.

The exported GenBank records include contig-level `source` features and `CDS` features for each coding gene call. When amino acid sequences are available in the contigs database, anvi'o also writes a `translation` qualifier for the CDS feature.

For example, to export only gene calls from a specific caller:

{{ codestart }}
anvi-export-antismash-input -c %(contigs-db)s \
                            -o output.gbk \
                            --gene-caller Prodigal
{{ codestop }}

If you want to skip short contigs while exporting:

{{ codestart }}
anvi-export-antismash-input -c %(contigs-db)s \
                            -o output.gbk \
                            --min-contig-length 5000
{{ codestop }}
