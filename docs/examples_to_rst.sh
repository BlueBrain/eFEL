# to call from root

# clean up
rm docs/source/sonata-network.rst
rm -rf docs/source/sonata-network_files
rm docs/source/nmc-portal.rst
rm -rf docs/source/L5TTPC2_files

# convert
jupyter nbconvert --to rst examples/sonata-network/sonata-network.ipynb
jupyter nbconvert --to rst examples/nmc-portal/L5TTPC2.ipynb

# move
mv examples/sonata-network/sonata-network.rst docs/source/
mv examples/sonata-network/sonata-network_files docs/source/
mv examples/nmc-portal/L5TTPC2.rst docs/source/nmc-portal.rst
mv examples/nmc-portal/L5TTPC2_files docs/source/