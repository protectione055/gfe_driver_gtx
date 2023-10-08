LD_LIBRARY_PATH=/usr/local/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/zhou822/BwGraph_Mono_Simple_Iterator/build/
export LD_LIBRARY_PATH 
for a in 4 8 16 24 32 40 48 56 64
do
    numactl -N 0 -l /home/zhou822/gfe_driver_bigdata/build/gfe_driver -G /home/zhou822/gfe_experiment_data/graph_source/graph500-22.properties -u -l bwgraph_rw -w $a
done