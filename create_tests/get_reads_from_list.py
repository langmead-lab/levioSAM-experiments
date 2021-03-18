import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--files', nargs='+',
        help='Lists of read names'
    )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Output prefix.'
    )
    parser.add_argument(
        '-1', '--raw_fastq1',
        help='Original FASTQ-R1 file.'
    )
    parser.add_argument(
        '-2', '--raw_fastq2',
        help='Original FASTQ-R2 file.'
    )
    return parser.parse_args()

def get_reads_from_list(files):
    lists = [[] for fn in files]
    for i, fn in enumerate(files):
        with open(fn, 'r') as f:
            for line in f:
                lists[i].append(line.rstrip())
    for i, fn in enumerate(files):
        lists[i] = set(lists[i])
        if i == 0:
            reads_isec = lists[i]
        else:
            reads_isec = reads_isec.intersection(lists[i])
    print(reads_isec)
    print(len(reads_isec))
    return reads_isec


def read_fastq(fn):
    read_dict = {}
    with open(fn, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.rstrip()
            if line_num % 4 == 0:
                read_dict[line[1:]] = [line]
                name = line[1:]
            else:
                read_dict[name].append(line)
    return read_dict


def process_io_fastq(fq1, fq2, reads_isec, output_prefix):
    read_dict1 = read_fastq(fq1)
    read_dict2 = read_fastq(fq2)
    f1 = open(output_prefix + '_1.fq', 'w')
    f2 = open(output_prefix + '_2.rq', 'w')
    for r in reads_isec:
        read1 = read_dict1[r]
        read2 = read_dict2[r]
        for l in read1:
            f1.write(l + '\n')
        for l in read2:
            f2.write(l + '\n')

if __name__ == '__main__':
    args = parse_args()
    reads_isec = get_reads_from_list(args.files)
    process_io_fastq(args.raw_fastq1, args.raw_fastq2, reads_isec, args.output_prefix)
