#!/usr/bin/env python3

import sys
import copy
import json

def transform_all(template_file, config):
    for entry in config:
        with open(entry["outpath"], "w") as outconn:
            transform(template_file, entry["replacements"], outconn)

def transform(template_file, replace_items, outconn):
    outfile = copy.deepcopy(template_file)
    for token, replacement in replace_items.items():
        outfile = outfile.replace(token, replacement)
    outconn.write(outfile)

def main():
    template_file = sys.stdin.read()
    with open(sys.argv[1], "r") as inconn:
        config = json.load(inconn)
    transform_all(template_file, config)

if __name__ == "__main__":
    main()
