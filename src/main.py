

import argparse
import os

from riptide_sampling import get_flux_samples

parser_general=argparse.ArgumentParser()

parser_general.add_argument("-d", "--dose-level", help="enter the dose level. 'Control','Low','Middle','High'",type=str)
parser_general.add_argument("-mol", "--molecule", help="enter the name of the molecule. amiodarone or valproic acid",type=str)
parser_general.add_argument("-rep", "--replicates", help="number of replicates. By default, riptide take the min and max over all replicates.",type=int)


args = parser_general.parse_args()

if args.molecule and args.dose_level:
    print(f"the {args.molecule} molecule is choosen.")
    print(f"the {args.dose_level} dose level is choosen")
    print(f"the {args.replicates} replicates is choosen")

    get_flux_samples('data/microarray/annotated_data_uniq_high_sd_flagged.tsv','data/microarray/open_tggates_cel_file_attribute.csv',sacrific_period='24 hr',dose_level=args.dose_level,compound_name=args.molecule, replicates=args.replicates)