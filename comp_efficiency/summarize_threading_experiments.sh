#python summarize_time_log_as_tsv.py -p single_end/lft-scaling -r 5 -t 1 2 4 8 16 32 64 -o single_end/lft-scaling.summary.tsv
#python summarize_time_log_as_tsv.py -p paired_end/lft-scaling -r 5 -t 1 2 4 8 16 32 64 -o paired_end/lft-scaling.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p single_end/lft -r 5 -t 1 2 4 8 16 32 64 -o single_end/lft-scaling.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p paired_end/lft -r 5 -t 1 2 4 8 16 32 64 -o paired_end/lft-scaling.summary.tsv
cat single_end/lft-scaling.summary.tsv > scaling.summary.tsv
#tail -n +2 paired_end/lft-scaling.summary.tsv | head -n -1 >> scaling.summary.tsv
tail -n +2 paired_end/lft-scaling.summary.tsv >> scaling.summary.tsv
