#!/bin/bash

./combine_goods <ex10_goods.txt | \
tee combined_goods.txt | \
./tukey_s.py tukey_prefix
