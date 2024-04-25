jupyter nbconvert --to rst ../examples/sonata-network/sonata-network.ipynb
jupyter nbconvert --to rst ../examples/nmc-portal/L5TTPC2.ipynb
cp ../examples/sonata-network/sonata-network.rst source/
cp -r ../examples/sonata-network/sonata-network_files source/
cp ../examples/nmc-portal/L5TTPC2.rst source/nmc-portal.rst
cp -r ../examples/nmc-portal/L5TTPC2_files source/