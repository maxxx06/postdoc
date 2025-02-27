

import riptide_sampling
import argparse

parser_general=argparse.ArgumentParser()
parser_general.add_argument('-f', nargs='?', default="config.yaml", type=str, help="path of your config file. defaut: config.yaml")

args = parser_general.parse_args()
riptide_sampling.load_riptide(args.f)
# riptide_sampling.load_analysis(args.f)
