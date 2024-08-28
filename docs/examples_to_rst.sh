# to call from root

# clean up
rm docs/source/sonata-network.rst
rm -rf docs/source/sonata-network_files
rm docs/source/nmc-portal.rst
rm -rf docs/source/L5TTPC2_files
rm docs/source/load_nwb.rst
rm -rf docs/source/load_nwb_files
rm docs/source/extrafeats_example.rst
rm docs/source/multiprocessing_example.rst
rm docs/source/voltage_clamp.rst
rm -rf docs/source/voltage_clamp_files

# convert
jupyter nbconvert --to rst examples/sonata-network/sonata-network.ipynb
jupyter nbconvert --to rst examples/nmc-portal/L5TTPC2.ipynb
jupyter nbconvert --to rst examples/neo/load_nwb.ipynb
jupyter nbconvert --to rst examples/extracellular/extrafeats_example.ipynb
jupyter nbconvert --to rst examples/parallel/multiprocessing_example.ipynb
jupyter nbconvert --to rst examples/voltage_clamp/voltage_clamp.ipynb

# move
mv examples/sonata-network/sonata-network.rst docs/source/
mv examples/sonata-network/sonata-network_files docs/source/
mv examples/nmc-portal/L5TTPC2.rst docs/source/nmc-portal.rst
mv examples/nmc-portal/L5TTPC2_files docs/source/
mv examples/neo/load_nwb.rst docs/source/
mv examples/neo/load_nwb_files docs/source/
mv examples/extracellular/extrafeats_example.rst docs/source/
mv examples/parallel/multiprocessing_example.rst docs/source/
mv examples/voltage_clamp/voltage_clamp.rst docs/source/
mv examples/voltage_clamp/voltage_clamp_files docs/source/