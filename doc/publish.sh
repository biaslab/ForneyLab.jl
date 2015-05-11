#!/bin/bash

# Push the html documentation to spsbrats.github.io
CHECKOUT="/tmp/spsbrats.github.io.git"

rm -f -r $CHECKOUT
git clone git@github.com:spsbrats/spsbrats.github.io.git $CHECKOUT
DOCPATH=$PWD
cd $CHECKOUT
git rm -f -r ForneyLab/documentation/*
mkdir --parents $CHECKOUT/ForneyLab/documentation
cp -r $DOCPATH/_build/html/* $CHECKOUT/ForneyLab/documentation/
git add ForneyLab/documentation/*
git commit -m "ForneyLab documentation update"
git push origin master
cd $DOCPATH
rm -f -r $CHECKOUT
echo "All done."