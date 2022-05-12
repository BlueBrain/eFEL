#!/bin/bash

set -e
set -x

git config --global http.proxy http://bbpproxy.epfl.ch:80/
git config --global https.proxy http://bbpproxy.epfl.ch:80/

tox_args='-v --recreate'
tox_args="${tox_args} -e py3-test"

if [ "${os}" = "bb5" ]
then
	. /opt/rh/rh-python36/enable
fi

cd $WORKSPACE

#########
# Virtualenv
#########

if [ ! -d "${WORKSPACE}/env" ]; then
  python -m venv env
fi

. ${WORKSPACE}/env/bin/activate
pip install pip --upgrade
pip install tox --upgrade   


#####
# Tests
#####

cd  ${WORKSPACE}/eFEL

tox ${tox_args}
