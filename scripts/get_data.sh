#!/usr/bin/env bash 
mkdir -p card/{canonical,prevalence}
wget -P card/canonical https://card.mcmaster.ca/latest/dat


https://card.mcmaster.ca/download/0/broadstreet-v3.0.3.tar.gz
wget -P card/prevalence https://card.mcmaster.ca/download/6/prevalence-v3.0.4.tar.gz

cd card/canonical
tar xvf broadstreet-v3.0.3.tar.gz

cd ../prevalence
tar xvf prevalence-v3.0.4.tar.gz

cd ../..
