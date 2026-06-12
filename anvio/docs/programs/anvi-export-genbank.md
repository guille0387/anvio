This program exports an anvi'o %(contigs-db)s as an annotated %(genbank-file)s.

The GenBank file produced by this program includes `CDS` features from anvi'o gene calls. This makes the output useful for downstream tools that need nucleotide sequences together with gene coordinates and protein translations.

In the exported GenBank file, each anvi'o coding gene call becomes a `CDS` feature. The anvi'o `gene_callers_id` is preserved in the `locus_tag` and `protein_id` qualifiers using this format:

{{ codestart }}
/locus_tag="anvio_gene_26"
/protein_id="anvio_gene_26"
{{ codestop }}

In this example, `anvio_gene_26` corresponds to `gene_callers_id` `26` in the %(contigs-db)s. This is a bridge that makes it possible for downstream analyses to map genes in the GenBank file back to anvi'o genes.

To export all contigs and all coding gene calls:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o %(genbank-file)s
{{ codestop }}

You can limit the export to contigs above a minimum length:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o %(genbank-file)s \
                    --min-contig-length 10000
{{ codestop }}

If your %(contigs-db)s contains gene calls from multiple sources, you can restrict the exported `CDS` features to one or more gene callers:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o %(genbank-file)s \
                    --gene-caller pyrodigal-gv
{{ codestop }}

You can provide multiple gene callers as a comma-separated list without spaces:

{{ codestart }}
anvi-export-genbank -c %(contigs-db)s \
                    -o %(genbank-file)s \
                    --gene-caller pyrodigal-gv,prodigal
{{ codestop }}

The resulting GenBank file can be used by external tools that accept annotated GenBank input. This program only exports the file; it does not run or import results from any downstream tool.
