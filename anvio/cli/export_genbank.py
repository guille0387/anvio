#!/usr/bin/env python

import sys

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

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['mateumenendez']
__requires__ = ['contigs-db']
__provides__ = ['genbank-file']
__description__ = "Export an anvi'o contigs database as an annotated GenBank file"


P = terminal.pluralize
pp = terminal.pretty_print


def sanity_check(args):
    """Make sure input and output paths are usable."""

    if not args.output_file:
        raise ConfigError("Please provide an output file path with `--output-file` so anvi'o knows where to write the "
                          "GenBank file.")

    filesnpaths.is_output_file_writable(args.output_file)
    filesnpaths.is_file_exists(args.contigs_db)

    if args.min_contig_length < 0:
        raise ConfigError("The minimum contig length must be 0 or greater.")


def get_gene_callers_id_token(gene_callers_id):
    """Return the stable token we put into GenBank qualifiers."""

    return f"anvio_gene_{gene_callers_id}"


def get_gene_calls_for_contig(contigs_db, gene_caller):
    """Return coding gene calls grouped by contig."""

    gene_calls_for_contig = {}
    gene_callers = gene_caller.split(',') if gene_caller else None
    available_gene_callers = sorted(set([gene_call['source'] for gene_call in contigs_db.genes_in_contigs_dict.values() if gene_call['call_type'] == constants.gene_call_types['CODING']]))

    if gene_callers:
        missing_gene_callers = sorted(set(gene_callers).difference(available_gene_callers))
        if missing_gene_callers:
            raise ConfigError(f"Anvi'o could not find {P('requested coding gene caller', len(missing_gene_callers))} "
                              f"in the contigs database: {', '.join(missing_gene_callers)}. Available coding gene "
                              f"callers are: {', '.join(available_gene_callers)}.")

    for gene_callers_id, gene_call in contigs_db.genes_in_contigs_dict.items():
        if gene_call['call_type'] != constants.gene_call_types['CODING']:
            continue

        if gene_callers and gene_call['source'] not in gene_callers:
            continue

        if gene_call['contig'] not in contigs_db.contig_sequences:
            continue

        contig = gene_call['contig']
        gene_calls_for_contig.setdefault(contig, [])
        gene_calls_for_contig[contig].append((gene_callers_id, gene_call))

    for contig in gene_calls_for_contig:
        gene_calls_for_contig[contig] = sorted(gene_calls_for_contig[contig], key=lambda item: (item[1]['start'], item[1]['stop']))

    return gene_calls_for_contig


def get_cds_feature_for_gene_call(gene_callers_id, gene_call, aa_sequences_dict):
    """Return a Biopython CDS feature for an anvi'o gene call."""

    strand = -1 if gene_call['direction'] == 'r' else 1
    location = FeatureLocation(int(gene_call['start']), int(gene_call['stop']), strand=strand)
    gene_callers_id_token = get_gene_callers_id_token(gene_callers_id)

    qualifiers = {
        'locus_tag': [gene_callers_id_token],
        'protein_id': [gene_callers_id_token],
        'codon_start': ['1'],
        'product': [f"anvi'o gene call {gene_callers_id}"],
        'note': [f"gene_callers_id={gene_callers_id};source={gene_call['source']}"],
    }

    if gene_call['partial']:
        qualifiers['note'][0] += ';partial=True'

    if gene_callers_id in aa_sequences_dict and aa_sequences_dict[gene_callers_id]['sequence']:
        qualifiers['translation'] = [aa_sequences_dict[gene_callers_id]['sequence']]

    return SeqFeature(location=location, type='CDS', qualifiers=qualifiers)


def get_source_feature_for_contig(sequence):
    """Return a source feature spanning the entire contig."""

    return SeqFeature(location=FeatureLocation(0, len(sequence), strand=1),
                      type='source',
                      qualifiers={'organism': ['unknown'],
                                  'mol_type': ['genomic DNA']})


def get_genbank_record_for_contig(contig_name, contig_sequence, gene_calls, aa_sequences_dict):
    """Return a Biopython SeqRecord for one contig."""

    record = SeqRecord(Seq(contig_sequence),
                       id=contig_name,
                       name=contig_name[:16],
                       description=f"{contig_name} exported from anvi'o")

    record.annotations['molecule_type'] = 'DNA'
    record.annotations['topology'] = 'linear'
    record.annotations['data_file_division'] = 'BCT'
    record.annotations['source'] = 'anvi-o'
    record.annotations['organism'] = 'unknown'

    record.features.append(get_source_feature_for_contig(contig_sequence))
    for gene_callers_id, gene_call in gene_calls:
        record.features.append(get_cds_feature_for_gene_call(gene_callers_id, gene_call, aa_sequences_dict))

    return record


def export_genbank(args, run, progress):
    """Write an annotated GenBank file."""

    sanity_check(args)

    quiet_run = terminal.Run(verbose=False)
    quiet_progress = terminal.Progress(verbose=False)

    contigs_db = dbops.ContigsSuperclass(args, r=quiet_run, p=quiet_progress)
    contigs_db.init_contig_sequences(min_contig_length=args.min_contig_length)

    db = dbops.ContigsDatabase(args.contigs_db)
    aa_sequences_dict = db.db.get_table_as_dict(t.gene_amino_acid_sequences_table_name)
    db.disconnect()

    if not contigs_db.genes_in_contigs_dict:
        raise ConfigError("This contigs database does not seem to have any gene calls. This exporter is specifically "
                          "for annotated GenBank output, so it needs gene calls from anvi'o.")

    gene_calls_for_contig = get_gene_calls_for_contig(contigs_db, args.gene_caller)
    if not gene_calls_for_contig:
        raise ConfigError("After filtering, anvi'o did not find any coding gene calls to report as CDS features. "
                          "Please check your contigs database, your `--gene-caller` filter, or your minimum contig "
                          "length.")

    records = []
    num_cds_features = 0
    num_cds_features_without_translation = 0

    progress.new('Preparing GenBank records', progress_total_items=len(contigs_db.contig_sequences))
    for contig_name in sorted(contigs_db.contig_sequences):
        progress.update(contig_name, increment=True)

        contig_sequence = contigs_db.contig_sequences[contig_name]['sequence']
        gene_calls = gene_calls_for_contig.get(contig_name, [])

        for gene_callers_id, _ in gene_calls:
            if gene_callers_id not in aa_sequences_dict or not aa_sequences_dict[gene_callers_id]['sequence']:
                num_cds_features_without_translation += 1

        num_cds_features += len(gene_calls)
        records.append(get_genbank_record_for_contig(contig_name, contig_sequence, gene_calls, aa_sequences_dict))

    progress.end()

    SeqIO.write(records, args.output_file, 'genbank')

    run.info('Contigs DB', args.contigs_db)
    run.info('Output GenBank', args.output_file)
    run.info('Contigs reported', pp(len(records)))
    run.info('CDS features reported', pp(num_cds_features))

    if num_cds_features_without_translation:
        run.warning(f"Anvi'o reported {P('CDS feature', num_cds_features_without_translation)} without a "
                    f"`translation` qualifier because amino acid sequences were not available in the contigs "
                    f"database for those gene calls. Some downstream tools may be unhappy about this.")


@time_program
def main():
    args = get_args()
    run = terminal.Run()
    progress = terminal.Progress()

    try:
        export_genbank(args, run, progress)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATABASE', "The anvi'o contigs database that will be translated into an annotated GenBank file.")
    groupA.add_argument(*anvio.A('contigs-db'), **anvio.K('contigs-db'))

    groupB = parser.add_argument_group('OUTPUT', "Where anvi'o should write the GenBank file.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))

    groupC = parser.add_argument_group('GENE CALL FILTERS', "Optional filters for gene calls and contigs.")
    groupC.add_argument(*anvio.A('gene-caller'), **anvio.K('gene-caller', {'help': "Which gene caller(s) would you like to include as CDS features? If providing multiple they should be comma-separated (no spaces). By default, all coding gene calls are included.", 'default': None,}))
    groupC.add_argument('-M', '--min-contig-length', default=0, type=int, metavar='INT',
                        help="Minimum contig length to include in the GenBank output. The default is %(default)d, "
                             "which reports every contig in the contigs database.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
