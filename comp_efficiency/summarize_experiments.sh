# python summarize_time_log_as_tsv.py -p single_end/bt2-t16 -r 5 -o single_end/bt2.summary.tsv
# python summarize_time_log_as_tsv.py -p single_end/lft-t16 -r 5 -o single_end/lft.summary.tsv
# python summarize_time_log_as_tsv.py -p paired_end/bt2-t16 -r 5 -o paired_end/bt2.summary.tsv
# python summarize_time_log_as_tsv.py -p paired_end/lft-t16 -r 5 -o paired_end/lft.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p single_end/bt2 -r 5 -o single_end/bt2.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p single_end/lft -r 5 -o single_end/lft.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p paired_end/bt2 -r 5 -o paired_end/bt2.summary.tsv
python ../scripts/summarize_time_log_as_tsv.py -p paired_end/lft -r 5 -o paired_end/lft.summary.tsv

cat single_end/bt2.summary.tsv > summary.tsv
tail -n +2 single_end/lft.summary.tsv | head -n -1 >> summary.tsv
tail -n +2 paired_end/bt2.summary.tsv | head -n -1 >> summary.tsv
tail -n +2 paired_end/lft.summary.tsv | head -n -1 >> summary.tsv
