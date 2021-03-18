import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input',
        help='Input raw SAM file'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output SAM file containing the target CIGAR operator'
    )
    parser.add_argument(
        '-p', '--op', default='S',
        help='Target CIGAR operator ["S"]'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    f_out = open(args.output, 'w')

    with open(args.input, 'r') as f:
        for line in f:
            if line[0] == '@':
                f_out.write(line)
            else:
                if line.split()[5].count(args.op) > 0:
                    f_out.write(line)
