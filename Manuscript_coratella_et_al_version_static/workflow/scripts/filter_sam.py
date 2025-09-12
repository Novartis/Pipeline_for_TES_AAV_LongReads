import sys, os, argparse



def primary_aln(sam_cols):
    # Return True if the alignment is the primary (bitwise flag = 0)
    # or the primary reverse complementary (bitwise flag = 16):
    if sam_cols[1] != '0' and sam_cols[1] != '16':
        return(True)
    else:
        return(False)


def bitflag_found(sam_cols, flags):
    # Return True if the bitwise flag is found in the inputted flags list.
    if int(sam_cols[1]) in flags:
        return(True)
    else:
        return(False)


def chr_chimera(sam_cols, sam_tags):
    if "SA" in sam_tags:
        for suppl_aln in sam_tags['SA'].split(';'):
            if suppl_aln == '':
                continue

            # If the 'SA' tag exists, get the chromosome of the supplementary alignment:
            sa_chr = suppl_aln.split(',')[0]

            # If the primary and supplementary alignments are not on the same chromosome, skip the read
            if sa_chr != sam_cols[2]:
                return(True)
    return(False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of SAM lines and applied various filters.""")

    # Input options:
    parser.add_argument("--chr_chimera", action='store_true',
        help="Filter chromosome chimeras by requiring primary and secondary alignments to be from the same chromosome.")

    parser.add_argument("--primary_aln", action='store_true',
        help="Filter alignments that are not primary alignments, defined as bitwise flag of 0 or 16.")

    parser.add_argument("--bitflag_excl", default=None, type=int, nargs='*',
        help="Filter out if any of the inputted bitflags are set.")

    parser.add_argument("--bitflag_incl", default=None, type=int, nargs='*',
        help="Only keep if any of the inputted bitflags are set.")

    parser.add_argument("--verbose", action='store_true')

    # Parse the arguments:
    args = parser.parse_args()

    if args.verbose:
        print(args, file=sys.stderr)


    # Read each SAM formatted line streamed from stdin:
    for sam_line in sys.stdin:
        # Skip header lines:
        sam_line = sam_line.strip()
        if sam_line.startswith('@'):
            print(sam_line)
            continue

        sam_cols = sam_line.split('\t')


        ### Apply SAM column filters ###
        if args.bitflag_excl:
            if bitflag_found(sam_cols, args.bitflag_excl):
                continue

        if args.bitflag_incl:
            if not bitflag_found(sam_cols, args.bitflag_incl):
                continue

        if args.primary_aln:
            if not primary_aln(sam_cols):
                continue


        # Read SAM tags if getting this far:
        sam_tags = dict()
        for tag in sam_cols[11:]:
            tag_split = tag.split(':')
            tag_val = ':'.join(tag_split[2:])
            if tag_split[1] == 'f':
                tag_val = float(tag_val)
            elif tag_split[1] == 'i':
                tag_val = int(tag_val)
            sam_tags[tag_split[0]] = tag_val


        ### Apply SAM tag filters ###
        if args.chr_chimera:
            if chr_chimera(sam_cols, sam_tags):
                continue


        # Getting past all the filters, print the line:
        print(sam_line)

