"""AmyloBeads: a tool for preparation and running of coarse-grained molecular
   dynamics of amyloidogenic polypeptides.

   (C) Egor Vasilenko, Dmitry Ivankov, Protein bioinformatics and evolution lab,
       Moscow Center for Molecular and Cellular Biology, 2025-2026.
"""


import argparse


parser = argparse.ArgumentParser(prog='AmyloBeads',
                                 description='''Prepare and run coarse-grained
                                                molecular dynamics of
                                                amyloidogenic polypeptides.''')

subparsers = parser.add_subparsers()

parser_extr = subparsers.add_parser('extr', help='''Extract geometry from
                                                    fibril PDB.''')
parser_extr.add_argument('--pdb', action='store_const',
                         help='PDB or mmCIF file with amyloid structure.')

parser_prep = subparsers.add_parser('prep', help='Prepare force field and cell.')

args = parser.parse_args()
print(args)
