import argparse
import re


SAM_NAME = 0
SAM_FLAG = 1
SAM_CHROM = 2
SAM_POS = 3
SAM_MAPQ = 4
SAM_CIGAR = 5
SAM_MCHROM = 6
SAM_MPOS = 7
SAM_TLEN = 8
SAM_SEQ = 9
SAM_QUAL = 10

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--files', nargs='+',
        help='Complete input SAM file (should be two files)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file. Supported ext: ".list", ".fastq"/".fq"'
    )
    parser.add_argument(
        '-r', '--raw_fastq',
        help='Original FASTQ file. Required if output is "*.fastq" or "*.fq"'
    )
    return parser.parse_args()


""" Returns a dict for a SAM file.

Output: a dict
    key: read name
    value: a list (should be of length=2) of lists (SAM fields from FLAG to QUAL).
"""
def read_sam(sam_fn):
    alns = {}
    with open(sam_fn, 'r') as f:
        for line in f:
            if line[0] != '@':
                line = line.split()
                alns.setdefault(line[0], []).append(line[: SAM_QUAL])
    return alns


def select_reads_core(list_alns, classify=True, max_class_count=30):
    invalid_names = []
    if classify:
        del_names = []
        ins_names = []
        start_s_names = []
        end_s_names = []
        other_names = []
    else:
        valid_names = []

    for name, fields in list_alns[0].items():
        if int(fields[0][SAM_FLAG]) & 64:
            first = fields[0]
            second = fields[1]
        else:
            first = fields[1]
            second = fields[0]
        queried_fields = list_alns[1][name]
        if int(queried_fields[0][SAM_FLAG]) & 64:
            q_first = queried_fields[0]
            q_second = queried_fields[1]
        else:
            q_first = queried_fields[1]
            q_second = queried_fields[0]
        
        if first == q_first and second == q_second:
            if not classify:
                valid_names.append(name)
            else:
                cigars = [first[SAM_CIGAR], second[SAM_CIGAR],
                          q_first[SAM_CIGAR], q_first[SAM_CIGAR]]
                # print(cigars)
                re_AtoZ = re.compile('[A-Z*]+')
                if max([c.count('D') for c in cigars]) > 0 and len(del_names) < max_class_count:
                    del_names.append(name)
                elif max([c.count('I') for c in cigars]) > 0 and len(ins_names) < max_class_count:
                    ins_names.append(name)
                elif sum([re_AtoZ.findall(c)[0] == 'S' for c in cigars]) > 0 and \
                    len(start_s_names) < max_class_count:
                    start_s_names.append(name)
                elif sum([re_AtoZ.findall(c)[-1] == 'S' for c in cigars]) > 0 and \
                    len(end_s_names) < max_class_count:
                    end_s_names.append(name)
                elif len(other_names) < max_class_count:
                    other_names.append(name)
        else:
            invalid_names.append(name)
            # print(name)
            print(first)
            print(second)
            print(q_first)
            print(q_second)
    print(len(invalid_names))
    if classify:
        print(len(del_names))
        print(len(ins_names))
        print(len(start_s_names))
        print(len(end_s_names))
        valid_names = del_names + ins_names + start_s_names + end_s_names + other_names
    else:
        print(len(valid_names))
    return valid_names


def select_reads_and_write(files, output, raw_fastq, output_ext):
    list_alns = [read_sam(fn) for fn in files]
    valid_names = select_reads_core(list_alns)
    if output_ext == 'list':
        with open(output, 'w') as f:
            for n in valid_names:
                f.write(n + '\n')
    print(valid_names)


if __name__ == '__main__':
    args = parse_args()
    output_ext = args.output.split('.')[-1]
    if output_ext == 'list':
        output_type = 'list'
    elif output_ext in ['fastq', 'fq']:
        output_type = 'fastq'
        assert args.raw_fastq
    assert len(args.files) == 2

    select_reads_and_write(
        args.files,
        args.output,
        args.raw_fastq,
        output_ext)
