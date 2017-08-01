#!/bin/sh

tests=`ls *.x`
for test in $tests; do
    echo "Running $test ..."
    ./$test
done
