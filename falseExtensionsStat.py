#!/usr/bin/env python
import numpy as np
import pandas as pd
import re
import json
import sys

inputFile = sys.argv[1]

rx="clones(?P<clones>[0-9]+)_coverage(?P<coverage>[0-9]+)_length(?P<length>[0-9]+)_seq(?P<_sequencer>[^_]+)"
results = pd.read_json(inputFile, lines=True)
columns_to_convert=["clones", "coverage", "length"]
results = pd.concat([results, results.inputFile.str.extract(rx, expand=True)], axis=1)
results[columns_to_convert] = results[columns_to_convert].apply(pd.to_numeric)
print("Total extension number: %s" % ((results.totalLExtensions.sum() + results.totalRExtensions.sum())))
print("False-extension number: %s" % (results.falseLExtensions.sum() + results.falseRExtensions.sum()))
print("False-extension percent total: %.4f%%" % (100.0 * (results.falseLExtensions.sum() + results.falseRExtensions.sum()) / (results.totalLExtensions.sum() + results.totalRExtensions.sum())))