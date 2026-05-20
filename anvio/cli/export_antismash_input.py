#!/usr/bin/env python

import sys

import anvio
import anvio.antismash as antismash

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import time_program


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['meren']
__requires__ = ['contigs-db']
__provides__ = ['genbank-file']
__description__ = "Export an anvi'o contigs database as annotated GenBank input for antiSMASH"


@time_program
def main():
    args = get_args()

    try:
        exporter = antismash.AntiSMASHInputExporter(args)
        exporter.process()
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-2)


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT DATABASE', "The anvi'o contigs database that will be translated into an antiSMASH-ready GenBank file.")
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
