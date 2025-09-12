"""
Copyright 2025 Novartis Institutes for BioMedical Research Inc.
 
Licensed under the MIT License (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
 
https://www.mit.edu/~amini/LICENSE.md
 
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import sys, os, argparse
import bisect



def primary_aln(sam_cols):
    # Return True if the alignment is the primary (bitwise flag = 0)
    # or the primary reverse complementary (bitwise flag = 16):
    if sam_cols[1] == '0' or sam_cols[1] == '16':
        return(True)
    else:
        return(False)


def bitflag_found(sam_cols, bitflags):
    # Return True if the bitwise flag is found in the inputted flags list.
    for flag in bitflags:
        if int(sam_cols[1]) & flag:
            return(True)
    return(False)


def bitflag_combi_found(sam_cols, flags):
    # Return True if the combined bitwise flag is found in the inputted flags list.
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


def secondary_chr(sam_cols, chr_name):
    if sam_cols[2] == chr_name:
        return(True)
    else:
        return(False)


def keep_Nseg(sam_line, sam_tags, cmd_args, sam_line_sorted):
    score = sam_tags[cmd_args[1]]
    # For descending order, we invert the value for sorting:
    if cmd_args[2].lower() == "dsc":
        score = -score

    # Add the line and trim:
    bisect.insort(sam_line_sorted, (score, sam_line))
    if len(sam_line_sorted) > cmd_args[0]:
        sam_line_sorted.pop()
    return(sam_line_sorted)


def keep_Nseg_chr(sam_line, sam_cols, sam_tags, cmd_args, sam_line_sorted):
    highest_aln_score = 1e6

    score = sam_tags[cmd_args[1]]
    assert(score/2 < highest_aln_score)

    # If the chromosome is matched keep the alignment:
    if sam_cols[2] == cmd_args[3]:
        score += highest_aln_score

    # For descending order, we invert the value for sorting:
    if cmd_args[2].lower() == "dsc":
        score = -score

    # Add the line and trim:
    bisect.insort(sam_line_sorted, (score, sam_line))
    if len(sam_line_sorted) > cmd_args[0]:
        sam_line_sorted.pop()
    return(sam_line_sorted)


def keep_Nseg_chr_chr(sam_line, sam_cols, sam_tags, cmd_args, sam_line_sorted):
    score = sam_tags[cmd_args[1]]
    # For descending order, we invert the value for sorting:
    if cmd_args[2].lower() == "dsc":
        score = -score

    # If the chromosome is matched keep the alignment:
    if sam_cols[2] == cmd_args[3]:
        # Add the line and trim:
        bisect.insort(sam_line_sorted[1], (score, sam_line))
        if len(sam_line_sorted[1]) > cmd_args[0]:
            sam_line_sorted[1].pop()
    else:
        # Add the line and trim:
        bisect.insort(sam_line_sorted[0], (score, sam_line))
        if len(sam_line_sorted[0]) > cmd_args[0]:
            sam_line_sorted[0].pop()
    return(sam_line_sorted)


def read_sam_tags(sam_cols):
    sam_tags = dict()
    for tag in sam_cols[11:]:
        tag_split = tag.split(':')
        tag_val = ':'.join(tag_split[2:])
        if tag_split[1] == 'f':
            tag_val = float(tag_val)
        elif tag_split[1] == 'i':
            tag_val = int(tag_val)
        sam_tags[tag_split[0]] = tag_val
    return(sam_tags)


def min_tag(sam_tags, tag, min_val):
    # Return True if the tag value is greater than or equal to the minimum value:
    if tag in sam_tags:
        if sam_tags[tag] >= min_val:
            return(True)
        else:
            return(False)
    else:
        return(False)

def max_tag(sam_tags, tag, max_val):
    # Return True if the tag value is less than or equal to the maximum value:
    if tag in sam_tags:
        if sam_tags[tag] <= max_val:
            return(True)
        else:
            return(False)
    else:
        return(False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""This takes a stream of SAM lines and applied various filters.""")

    # Input options:
    parser.add_argument("--chr_chimera", action='store_true',
        help="Filter chromosome chimeras by requiring primary and supplementary alignments to be from the same chromosome.")

    parser.add_argument("--secondary_chr", default=None, type=str,
        help="Filter out non-primary alignments if segment is _not_ aligned to specified chromosome name.")

    parser.add_argument("--primary_aln", action='store_true',
        help="Filter alignments that are not primary alignments, defined as bitwise flag of 0 or 16.")

    parser.add_argument("--bitflag_combi_excl", default=None, type=int, nargs='*',
        help="Filter out if any of the inputted bitflag combinations are set.")

    parser.add_argument("--bitflag_combi_incl", default=None, type=int, nargs='*',
        help="Only keep if any of the inputted bitflag combinations are set.")

    parser.add_argument("--bitflag_excl", default=None, type=int, nargs='*',
        help="Filter out if any of the inputted bitflag are set.")

    parser.add_argument("--bitflag_incl", default=None, type=int, nargs='*',
        help="Only keep if any of the inputted bitflag are set.")

    parser.add_argument("--min_tag", default=None, nargs='+', action='append', metavar=('MIN', 'TAG'),
        help="Only keep if tag is minimum this value. If tag not found also discard the line.")

    parser.add_argument("--max_tag", default=None, nargs='+', action='append', metavar=('MAX', 'TAG'),
        help="Only keep if tag is maximum this value. If tag not found also discard the line.")

    parser.add_argument("--keep_Nseg", default=None, nargs=3, metavar=('INTEGER', 'STRING', 'STRING'),
        help="""After filtering, keep the first N aligned segments of a query sequence
                after sorting according to a SAM tag, either ascending or descending.
                Specify as: interger TAG asc/dsc e.g. 2 AS dsc.
                Expects the SAM stream to be sorted by query name i.e. `samtools sort -n`.""")

    parser.add_argument("--keep_Nseg_chr", default=None, nargs=4, metavar=('INTEGER', 'STRING', 'STRING', 'STRING'),
        help="""Same as --keep_Nseg but with one extra argument for a chromosome
                to bump up to the top rank.
                Specify as: interger TAG asc/dsc chr e.g. 2 AS dsc AAV6_vector.""")

    parser.add_argument("--keep_Nseg_chr_chr", default=None, nargs=4, metavar=('INTEGER', 'STRING', 'STRING', 'STRING'),
        help="""Same as --keep_Nseg but sorting into two piles defined by the extra argument specifying a chromosome.
                One pile for segments mapped to the specified chromosome, another pile for the rest.
                Notice, that the first argument integer specifies how many segments to extract from each pile
                i.e. if 1 is specified, at most two segments are returned.
                Specify as: interger TAG asc/dsc chr e.g. 1 AS dsc AAV6_vector.""")

    parser.add_argument("--ref_incl", default=None, nargs='*', type=str,
        help="Only keep if the reference name is in the list. If not found, discard the line.")

    parser.add_argument("--ref_excl", default=None, nargs='*', type=str,
        help="Only keep if the reference name is not in the list. If found, discard the line.")

    parser.add_argument("--verbose", action='store_true')

    # Parse the arguments:
    args = parser.parse_args()

    # Segment filters are mutually exclusive:
    if sum(1 for b in [args.keep_Nseg, args.keep_Nseg_chr, args.keep_Nseg_chr_chr] if b) > 1:
        raise Exception("Command line arguments --keep_Nseg, --keep_Nseg_chr and --keep_Nseg_chr_chr are mutually exclusive.")
    elif args.keep_Nseg or args.keep_Nseg_chr:
        sam_line_sorted = list()
        last_qname = ''
    elif args.keep_Nseg_chr_chr:
        sam_line_sorted = [[], []]
        last_qname = ''

    # First argument for segment filters is an integer:
    for seg_filt in ["keep_Nseg", "keep_Nseg_chr", "keep_Nseg_chr_chr"]:
        if vars(args)[seg_filt]:
            try:
                vars(args)[seg_filt][0] = int(vars(args)[seg_filt][0])
            except ValueError:
                parser.error(f"The first value for --{seg_filt} must be an integer.")

    # Convert reference include/exclude lists to sets for faster lookup:
    if args.ref_incl:
        args.ref_incl = set(args.ref_incl)
    if args.ref_excl:
        args.ref_excl = set(args.ref_excl)

    if args.verbose:
        print(args, file=sys.stderr)


    
    # Read each SAM formatted line streamed from stdin:
    sam_line = sys.stdin.readline()
    while sam_line:
        # Skip header lines:
        sam_line = sam_line.strip()
        if sam_line.startswith('@'):
            print(sam_line)
            sam_line = sys.stdin.readline()
            continue

        sam_cols = sam_line.split('\t')


        ### Apply SAM column filters ###
        if args.bitflag_excl:
            if bitflag_found(sam_cols, args.bitflag_excl):
                sam_line = sys.stdin.readline()
                continue

        if args.bitflag_incl:
            if not bitflag_found(sam_cols, args.bitflag_incl):
                sam_line = sys.stdin.readline()
                continue

        if args.bitflag_combi_excl:
            if bitflag_combi_found(sam_cols, args.bitflag_combi_excl):
                sam_line = sys.stdin.readline()
                continue

        if args.bitflag_combi_incl:
            if not bitflag_combi_found(sam_cols, args.bitflag_combi_incl):
                sam_line = sys.stdin.readline()
                continue

        if args.primary_aln:
            if not primary_aln(sam_cols):
                sam_line = sys.stdin.readline()
                continue

        if args.secondary_chr and not primary_aln(sam_cols):
            if not secondary_chr(sam_cols, args.secondary_chr):
                sam_line = sys.stdin.readline()
                continue

        if args.ref_incl:
            if sam_cols[2] not in args.ref_incl:
                sam_line = sys.stdin.readline()
                continue

        if args.ref_excl:
            if sam_cols[2] in args.ref_excl:
                sam_line = sys.stdin.readline()
                continue


        # Read SAM tags if required:
        if (args.chr_chimera or args.keep_Nseg or args.keep_Nseg_chr or args.keep_Nseg_chr_chr
            or args.min_tag or args.max_tag):
            sam_tags = read_sam_tags(sam_cols)

        ### Apply SAM tag filters ###
        if args.chr_chimera:
            if chr_chimera(sam_cols, sam_tags):
                sam_line = sys.stdin.readline()
                continue
        
        if args.min_tag:
            passed = True
            for mta in args.min_tag:
                try:
                    min_val = float(mta[0])
                except ValueError:
                    parser.error(f"The first value for --min_tag must be a number: {mta[0]}")
                tag = mta[1]
                if not min_tag(sam_tags, tag, min_val):
                    passed = False
            if not passed:
                sam_line = sys.stdin.readline()
                continue

        if args.max_tag:
            passed = True
            for mta in args.max_tag:
                try:
                    max_val = float(mta[0])
                except ValueError:
                    parser.error(f"The first value for --max_tag must be a number: {mta[0]}")
                tag = mta[1]
                if not max_tag(sam_tags, tag, max_val):
                    passed = False
            if not passed:
                sam_line = sys.stdin.readline()
                continue


        ### Apply segment filters ###
        if args.keep_Nseg or args.keep_Nseg_chr or args.keep_Nseg_chr_chr:
            # New query:
            if sam_cols[0] != last_qname:
                # Update query name:
                last_qname = sam_cols[0]

                # Dump saved sam lines for the query.
                # If keeping two ordered lists:
                if args.keep_Nseg_chr_chr:
                    # Dump list one:
                    for slt in sam_line_sorted[0]:
                        print(slt[1])
                    # Dump list two:
                    for slt in sam_line_sorted[1]:
                        print(slt[1])
                    sam_line_sorted = [[], []]

                # If keep a single ordered list:
                else:
                    for slt in sam_line_sorted:
                        print(slt[1])
                    sam_line_sorted = list()

            # Put sam line into sorted list:
            if args.keep_Nseg:
                sam_line_sorted = keep_Nseg(sam_line, sam_tags, args.keep_Nseg, sam_line_sorted)
            elif args.keep_Nseg_chr:
                sam_line_sorted = keep_Nseg_chr(sam_line, sam_cols, sam_tags, args.keep_Nseg_chr, sam_line_sorted)
            elif args.keep_Nseg_chr_chr:
                sam_line_sorted = keep_Nseg_chr_chr(sam_line, sam_cols, sam_tags, args.keep_Nseg_chr_chr, sam_line_sorted)
        else:
            # No segment filters, passed other filters, print the line:
            print(sam_line)

        # Update sam line for next round:
        sam_line = sys.stdin.readline()

    # Dump saved sam lines for the last query.
    if args.keep_Nseg_chr_chr:
        for slt in sam_line_sorted[0]:
            print(slt[1])
        for slt in sam_line_sorted[1]:
            print(slt[1])
    elif args.keep_Nseg or args.keep_Nseg_chr:
        for slt in sam_line_sorted:
            print(slt[1])






