"""Summarize time_log files as a TSV.

Example: python summarize_time_log_as_tsv -p single_end/lft -r 5 -o single_end/lft.summary.tsv
"""
import argparse
import numpy as np
import pandas as pd

def file_parsing(
        fn, label, list_usr_time, list_sys_time, 
        list_percent_cpu, list_wall_time, list_max_res_size, list_label):
    with open(fn, 'r') as f:
        for line in f:
            if line.count('User time (seconds):'):
                list_usr_time.append(float(line.split(':')[1]))
            elif line.count('System time (seconds)'):
                list_sys_time.append(float(line.split(':')[1]))
            elif line.count('Percent of CPU this job got:'):
                list_percent_cpu.append(int(line.split(':')[1][:-2]))
            elif line.count('Elapsed (wall clock) time (h:mm:ss or m:ss):'):
                wt_raw = line.split('):')[1]
                wt_raw = [float(t) for t in wt_raw.split(':')]
                if len(wt_raw) == 3:
                    wt = wt_raw[0] * 3600 + wt_raw[1] * 60 + wt_raw[2]
                else:
                    wt = wt_raw[0] * 60 + wt_raw[1]
                list_wall_time.append(wt)
            elif line.count('Maximum resident set size (kbytes):'):
                list_max_res_size.append(int(line.split(':')[1]))
        list_label.append(label)

def summarize_time_log_as_tsv(prefix, num_reps, threads_exp, output):
    list_usr_time = []
    list_sys_time = []
    list_percent_cpu = []
    list_wall_time = []
    list_max_res_size = []
    list_label = []
    list_thread = []
    for i in range(1, 1 + num_reps):
        if len(threads_exp) == 1:
            fn = f'{prefix}-{i}.time_log'
            file_parsing(fn, prefix, list_usr_time, list_sys_time,
                    list_percent_cpu, list_wall_time, list_max_res_size, list_label)
            list_thread.append(threads_exp[0])
        else:
            for t in threads_exp:
                #fn = f'{prefix}-t{t}.time_log'
                fn = f'{prefix}-{i}-t{t}.time_log'
                file_parsing(fn, f'{prefix}', list_usr_time, list_sys_time,
                        list_percent_cpu, list_wall_time, list_max_res_size, list_label)
                list_thread.append(t)
    list_cpu_time = [u + list_sys_time[i] for i, u in enumerate(list_usr_time)]
    # cpu_time_sd = np.std(list_cpu_time)
    # percent_cpu_sd = np.std(list_percent_cpu)
    # wall_time_sd = np.std(list_wall_time)
    # max_res_size_sd = np.std(list_max_res_size)
    df = pd.DataFrame(
            {"CpuTime": list_cpu_time, 
                # "SD-CpuTime": cpu_time_sd,
                "UsrTime": list_usr_time, "SysTime": list_sys_time,
                "PrecentCpu": list_percent_cpu, 
                # "SD-PercentCpu": percent_cpu_sd,
                "WallTime": list_wall_time, 
                # "SD-WallTime": wall_time_sd,
                "MaxResSize": list_max_res_size, 
                # "SD-MaxResSize": max_res_size_sd,
                "Label": list_label,
                "Threads": list_thread})
    df.to_csv(output, sep='\t', index=None)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--prefix',
        help='Prefix of time_log files'
    )
    parser.add_argument(
        '-r', '--num_reps', type=int, default=5,
        help='Number of replica [5]'
    )
    parser.add_argument(
        '-t', '--threads_exp', nargs='+', type=int, default=[16],
        help='Thread experiment parameters. Example: -t 1 2 4 8. [16]'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file name'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    summarize_time_log_as_tsv(args.prefix, args.num_reps, args.threads_exp, args.output)
