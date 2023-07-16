#!/bin/sh
LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/zhou822/BwGraph_v2/build/
export LD_LIBRARY_PATH
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 64 -d bwgraph_yahoo_ordered_64.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 56 -d bwgraph_yahoo_ordered_56.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 48 -d bwgraph_yahoo_ordered_48.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 40 -d bwgraph_yahoo_ordered_40.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 32 -d bwgraph_yahoo_ordered_32.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 24 -d bwgraph_yahoo_ordered_24.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 16 -d bwgraph_yahoo_ordered_16.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 12 -d bwgraph_yahoo_ordered_12.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 8 -d bwgraph_yahoo_ordered_8.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 4 -d bwgraph_yahoo_ordered_4.sqlite3 --is_timestamped true