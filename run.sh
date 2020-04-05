#!/usr/bin/env bash

rm input.xy output.dat
python generate.py $1
python md.py $1

