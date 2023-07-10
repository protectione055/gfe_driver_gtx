#!/bin/sh
LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/zhou822/BwGraph_v2/build/
export LD_LIBRARY_PATH
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 64 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 56 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 40 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 32 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 24 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 16 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 12 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 8 -d output_results.sqlite3 --is_timestamped true
/home/zhou822/gfe_driver_sortledton/build/gfe_driver -G /home/zhou822/gfe_experiment_data/yahoo-song/out.yahoo-song.el -u -l bwgraph_rw -w 4 -d output_results.sqlite3 --is_timestamped true