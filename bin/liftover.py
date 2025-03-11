from pyliftover import LiftOver
import pandas as pd
import argparse
import sys

def main(map_file, chain_file, out_prefix):
    map = pd.read_csv(map_file, sep='\t', header=None, names=["CHR", "SNP", "CM", "BP"])
    lo = LiftOver(chain_file)
    unmapped = []
    new_map = []
    for i, row in map.iterrows():
        chr_str = row.CHR if str(row.CHR).startswith("chr") else f"chr{row.CHR}"
        result = lo.convert_coordinate(chr_str, int(row.BP))
        if result is None or len(result) != 1:
            new_map.append([row.CHR, row.SNP, row.CM, row.BP])
            unmapped.append(row.SNP)
        else:
            new_map.append([result[0][0], row.SNP, row.CM, result[0][1]])

    new_map = pd.DataFrame(new_map, columns=["CHR", "SNP", "CM", "BP"])
    # Write the new map file
    new_map.to_csv(f"{out_prefix}.map", sep='\t', header=False, index=False)
    # Write the unmapped SNPs
    with open(f"{out_prefix}.unmapped", 'w') as f:
        f.write("\n".join(unmapped))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="%(prog)s converts bim file to another, using the provided chain file."
    )
    parser.add_argument('-m', "--map", dest='mapFile', required = True,
                        help='The plink MAP file to liftOver.')
    parser.add_argument('-o', "--out", dest='outPrefix', required = True,
                        help='The prefix to give to the output files.')
    parser.add_argument('-c', "--chain", dest='chainFile', required = True,
                        help='The location of the chain file')

    # Show usage message if user hasn't provided any arguments, rather
    # than giving a non-descript error message with the usage()
    if len(sys.argv) == 1:
      parser.print_help()
      sys.exit()

    args = parser.parse_args()

    main(args.mapFile, args.chainFile, args.outPrefix)


