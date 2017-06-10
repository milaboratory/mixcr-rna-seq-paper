#!/usr/bin/env python
import numpy as np
import pandas as pd
import re
import json
import sys

readToDescrFile = sys.argv[1]
overlapped = sys.argv[2]
assembleReport = sys.argv[3]

readIdToCDR3 = pd.read_table(readToDescrFile)
readIdToCDR3["CDR3"] = readIdToCDR3.descrR1.str.extract("GClone\|([^|]*)\|", expand=False)
del readIdToCDR3["descrR1"]
rescued = pd.read_table(overlapped)
rescued["targetDescriptionsNoBrackets"] = rescued.targetDescriptions.str.replace(r'\[[^\[\]]*\]', '')
rx="VJOverlap\([0-9]+\) = [LR](?P<R1>[0-9]+)\.[01] \+ [LR](?P<R2>[0-9]+)\.[01]"
rescued = pd.concat([rescued, rescued.targetDescriptionsNoBrackets.str.extract(rx, expand=True).apply(pd.to_numeric)], axis=1)
m1 = rescued.merge(readIdToCDR3, left_on="R1", right_on="readId", how="inner")
m2 = m1.merge(readIdToCDR3, left_on="R2", right_on="readId", how="inner")

clones = 0
for line in open(assembleReport,'r'):
    m = re.search("Reads used in clonotypes, percent of total: ([0-9]+)", line)
    if m:
        clones = int(m.group(1))
        break

result = {
    "inputFile": overlapped,
    "clonesTotal": clones,
    "totalAlignments": len(rescued),
    "totalAlignmentsWithCDR3": int(rescued.nSeqCDR3.notnull().sum()),
    "totalOverlaps": len(m2),
    "hqOverlaps": int((m2.minQualCDR3.fillna(0.0) >= 20).sum()),
    "correctOverlaps": int((m2.CDR3_x == m2.CDR3_y).sum()),
    "overlapsFromDifferentClones": int((m2.CDR3_x != m2.CDR3_y).sum()),
    "overlapsProducingNewCDR3": int(((m2.CDR3_x != m2.nSeqCDR3) & (m2.CDR3_y != m2.nSeqCDR3) & (m2.CDR3_x != m2.CDR3_y)).sum()), 
    "hqOverlapsProducingNewCDR3": int(((m2.minQualCDR3.fillna(0.0) >= 20) & (m2.CDR3_x != m2.nSeqCDR3) & (m2.CDR3_y != m2.nSeqCDR3) & (m2.CDR3_x != m2.CDR3_y)).sum()),
    "newCDR3Diversity": len(m2.loc[(m2.CDR3_x != m2.nSeqCDR3) & (m2.CDR3_y != m2.nSeqCDR3) & (m2.CDR3_x != m2.CDR3_y), "nSeqCDR3"].unique()), 
    "hqNewCDR3Diversity": len(m2.loc[(m2.minQualCDR3.fillna(0.0) >= 20) & (m2.CDR3_x != m2.nSeqCDR3) & (m2.CDR3_y != m2.nSeqCDR3) & (m2.CDR3_x != m2.CDR3_y), "nSeqCDR3"].unique()) 
}
print(json.dumps(result))
