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
results["percentBad"] = results.overlapsProducingNewCDR3 / results.totalOverlaps
results["percentBadDiversity"] = results.newCDR3Diversity / results.clonesTotal
results["percentBadDiversityHQ"] = results.hqNewCDR3Diversity / results.clonesTotal
print("Bad overlaps in 10^2 and 10^3: %.3f%%" % (100.0 * results.loc[results.clones <= 1000, "percentBad"].max()))
print("Bad overlap diversity in 10^2 and 10^3: %.3f%%" % (100.0 * results.loc[results.clones <= 1000, "percentBadDiversity"].max()))
print("Bad overlap high quality diversity in 10^2 and 10^3: %.3f%%" % (100.0 * results.loc[results.clones <= 1000, "percentBadDiversityHQ"].max()))
print("Bad overlaps: %.3f%%" % (100.0 * results.percentBad.max()))
print("Bad overlaps diversity: %.3f%%" % (100.0 * results.percentBadDiversity.max()))
print("Bad overlaps high quality diversity: %.3f%%" % (100.0 * results.percentBadDiversityHQ.max()))