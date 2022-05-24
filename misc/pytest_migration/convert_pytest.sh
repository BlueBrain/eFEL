#!/bin/bash

cd efel/tests

declare -a StringArray=("./")

for dir in ${StringArray[@]}
do
    sed -i'' 's/nt.assert_raises/pytest.raises/g' $dir*.py
    sed -i'' 's/nt.ok_/nt.assert_true/g' $dir*.py
    sed -i'' 's/nt.eq_/nt.assert_equal/g' $dir*.py
    sed -i'' 's/nt.assert/assert/g' $dir*.py
    sed -i'' 's/nt.assert_almost_equal/numpy.testing.assert_allclose/g' $dir*.py
    sed -i'' 's/places=/rtol=0, atol=1e-/g' $dir*.py
    sed -i'' 's/@nt.raises(/@pytest.mark.xfail(raises=/g' $dir*.py
    sed -i'' 's/import nose.tools as nt//g' $dir*.py 
done

nose2pytest -v .