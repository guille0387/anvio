#!/usr/bin/env python
"""Utilities for working with antiSMASH from anvi'o."""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import anvio
import anvio.dbops as dbops
import anvio.tables as t
import anvio.constants as constants
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "The Anvi'o Project"
__email__ = "a.murat.eren@gmail.com"


run = terminal.Run()
progress = terminal.Progress()
pp = terminal.pretty_print
P = terminal.pluralize


class AntiSMASHInputExporter(object):
    """Generate antiSMASH-ready GenBank input from an anvi'o contigs database."""

    def __init__(self, args, run=run, progress=progress):
        self.args = args
        self.run = run
        self.progress = progress

        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.contigs_db_path = A('contigs_db')
        self.output_file_path = A('output_file')
        self.gene_caller = A('gene_caller')
        self.min_contig_length = A('min_contig_length') or 0

        self.contigs_db = None
        self.aa_sequences_dict = {}

        self.sanity_check()


    def sanity_check(self):
        """Make sure input and output paths are usable."""

        if not self.output_file_path:
            raise ConfigError("Please provide an output file path with `--output-file` so anvi'o knows where to "
                              "write the GenBank file for antiSMASH.")

        filesnpaths.is_output_file_writable(self.output_file_path)
        filesnpaths.is_file_exists(self.contigs_db_path)

        if self.min_contig_length < 0:
            raise ConfigError("The minimum contig length must be 0 or greater.")


    def init_contigs_db(self):
        """Initialize the contigs database and its sequence data."""

        quiet_run = terminal.Run(verbose=False)
        quiet_progress = terminal.Progress(verbose=False)

        self.contigs_db = dbops.ContigsSuperclass(self.args, r=quiet_run, p=quiet_progress)
        self.contigs_db.init_contig_sequences(min_contig_length=self.min_contig_length)

        db = dbops.ContigsDatabase(self.contigs_db_path)
        self.aa_sequences_dict = db.db.get_table_as_dict(t.gene_amino_acid_sequences_table_name)
        db.disconnect()


    def get_gene_callers_id_token(self, gene_callers_id):
        """Return the stable token we put into GenBank qualifiers."""

        return f"anvio_gene_{gene_callers_id}"


    def get_gene_calls_for_contig(self):
        """Return coding gene calls grouped by contig."""

        gene_calls_for_contig = {}
        gene_callers = self.gene_caller.split(',') if self.gene_caller else None
        available_gene_callers = sorted(set([gene_call['source'] for gene_call in self.contigs_db.genes_in_contigs_dict.values() if gene_call['call_type'] == constants.gene_call_types['CODING']]))

        if gene_callers:
            missing_gene_callers = sorted(set(gene_callers).difference(available_gene_callers))
            if missing_gene_callers:
                raise ConfigError(f"Anvi'o could not find {P('requested coding gene caller', len(missing_gene_callers))} "
                                  f"in the contigs database: {', '.join(missing_gene_callers)}. Available coding gene "
                                  f"callers are: {', '.join(available_gene_callers)}.")

        for gene_callers_id, gene_call in self.contigs_db.genes_in_contigs_dict.items():
            if gene_call['call_type'] != constants.gene_call_types['CODING']:
                continue

            if gene_callers and gene_call['source'] not in gene_callers:
                continue

            if gene_call['contig'] not in self.contigs_db.contig_sequences:
                continue

            contig = gene_call['contig']
            gene_calls_for_contig.setdefault(contig, [])
            gene_calls_for_contig[contig].append((gene_callers_id, gene_call))

        for contig in gene_calls_for_contig:
            gene_calls_for_contig[contig] = sorted(gene_calls_for_contig[contig], key=lambda item: (item[1]['start'], item[1]['stop']))

        return gene_calls_for_contig


    def get_cds_feature_for_gene_call(self, gene_callers_id, gene_call):
        """Return a Biopython CDS feature for an anvi'o gene call."""

        strand = -1 if gene_call['direction'] == 'r' else 1
        location = FeatureLocation(int(gene_call['start']), int(gene_call['stop']), strand=strand)
        gene_callers_id_token = self.get_gene_callers_id_token(gene_callers_id)

        qualifiers = {
            'locus_tag': [gene_callers_id_token],
            'protein_id': [gene_callers_id_token],
            'codon_start': ['1'],
            'product': [f"anvi'o gene call {gene_callers_id}"],
            'note': [f"gene_callers_id={gene_callers_id};source={gene_call['source']}"],
        }

        if gene_call['partial']:
            qualifiers['note'][0] += ';partial=True'

        if gene_callers_id in self.aa_sequences_dict and self.aa_sequences_dict[gene_callers_id]['sequence']:
            qualifiers['translation'] = [self.aa_sequences_dict[gene_callers_id]['sequence']]

        return SeqFeature(location=location, type='CDS', qualifiers=qualifiers)


    def get_source_feature_for_contig(self, sequence):
        """Return a source feature spanning the entire contig."""

        return SeqFeature(location=FeatureLocation(0, len(sequence), strand=1),
                          type='source',
                          qualifiers={'organism': ['unknown'],
                                      'mol_type': ['genomic DNA']})


    def get_genbank_record_for_contig(self, contig_name, contig_sequence, gene_calls):
        """Return a Biopython SeqRecord for one contig."""

        record = SeqRecord(Seq(contig_sequence),
                           id=contig_name,
                           name=contig_name[:16],
                           description=f"{contig_name} exported from anvi'o for antiSMASH")

        record.annotations['molecule_type'] = 'DNA'
        record.annotations['topology'] = 'linear'
        record.annotations['data_file_division'] = 'BCT'
        record.annotations['source'] = 'anvi-o'
        record.annotations['organism'] = 'unknown'

        record.features.append(self.get_source_feature_for_contig(contig_sequence))
        for gene_callers_id, gene_call in gene_calls:
            record.features.append(self.get_cds_feature_for_gene_call(gene_callers_id, gene_call))

        return record


    def process(self):
        """Write an antiSMASH-compatible GenBank file."""

        self.init_contigs_db()

        if not self.contigs_db.genes_in_contigs_dict:
            raise ConfigError("This contigs database does not seem to have any gene calls. antiSMASH can run on "
                              "plain FASTA, but this exporter is specifically for annotated GenBank input, so it "
                              "needs gene calls from anvi'o.")

        gene_calls_for_contig = self.get_gene_calls_for_contig()
        if not gene_calls_for_contig:
            raise ConfigError("After filtering, anvi'o did not find any coding gene calls to report as CDS features. "
                              "Please check your contigs database, your `--gene-caller` filter, or your minimum contig "
                              "length.")

        records = []
        num_cds_features = 0
        num_cds_features_without_translation = 0

        self.progress.new('Preparing GenBank records', progress_total_items=len(self.contigs_db.contig_sequences))
        for contig_name in sorted(self.contigs_db.contig_sequences):
            self.progress.update(contig_name, increment=True)

            contig_sequence = self.contigs_db.contig_sequences[contig_name]['sequence']
            gene_calls = gene_calls_for_contig.get(contig_name, [])

            for gene_callers_id, _ in gene_calls:
                if gene_callers_id not in self.aa_sequences_dict or not self.aa_sequences_dict[gene_callers_id]['sequence']:
                    num_cds_features_without_translation += 1

            num_cds_features += len(gene_calls)
            records.append(self.get_genbank_record_for_contig(contig_name, contig_sequence, gene_calls))

        self.progress.end()

        SeqIO.write(records, self.output_file_path, 'genbank')

        self.run.info('Contigs DB', self.contigs_db_path)
        self.run.info('Output GenBank', self.output_file_path)
        self.run.info('Contigs reported', pp(len(records)))
        self.run.info('CDS features reported', pp(num_cds_features))

        if num_cds_features_without_translation:
            self.run.warning(f"Anvi'o reported {P('CDS feature', num_cds_features_without_translation)} without a "
                             f"`translation` qualifier because amino acid sequences were not available in the contigs "
                             f"database for those gene calls. If antiSMASH complains about missing translations, this "
                             f"will be the first place to look.")
