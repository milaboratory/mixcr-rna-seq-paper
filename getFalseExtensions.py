#!/usr/bin/env python
import numpy as np
import pandas as pd
import re
import json
import sys

inputFile = sys.argv[1]
assembleReport = sys.argv[2]

extensions = pd.read_table(inputFile)
numericColumns = ["refPoints","readSequence", "targetDescriptions"]
otherColumns = ["allVHitsWithScore","allJHitsWithScore"]
extensionsR1 = extensions.descrR1.rename("descr").to_frame()

for c in numericColumns:
    extensionsR1[c] = extensions[c].str.split(",").apply(lambda d: d[0])
for c in otherColumns:
    extensionsR1[c] = extensions[c]

extensionsR2 = extensions.descrR1.rename("descr").to_frame()
for c in numericColumns:
    extensionsR2[c] = extensions[c].str.split(",").apply(lambda d: d[1] if len(d) > 1 else "")
for c in otherColumns:
    extensionsR2[c] = extensions[c]

extensionsR2 = extensionsR2[extensionsR2.refPoints != ""]
extensionsFull = extensionsR1.append(extensionsR2, ignore_index=True)
anchorPointsRegex="^(?:-?[0-9]*:){8}(?:-?[0-9]*):(?P<CDR3Begin>-?[0-9]*):(?P<V3Deletion>-?[0-9]*):(?P<VEnd>-?[0-9]*):(?P<DBegin>-?[0-9]*):(?P<D5Deletion>-?[0-9]*):(?P<D3Deletion>-?[0-9]*):(?P<DEnd>-?[0-9]*):(?P<JBegin>-?[0-9]*):(?P<J5Deletion>-?[0-9]*):(?P<CDR3End>-?[0-9]*):(?:-?[0-9]*:){2}(?:-?[0-9]*)$"
extensionsFull = pd.concat([extensionsFull, extensionsFull.refPoints.str.extract(anchorPointsRegex, expand=True).apply(pd.to_numeric)], axis=1)
extensionsFull["trueCDR3"] = extensionsFull.descr.str.split("|").apply(lambda d: d[1])
extensionsFull["topVScore"] = extensionsFull.allVHitsWithScore.str.extract("^[^(]+\((?P<extended>[0-9.]+)\)", expand=False).apply(pd.to_numeric)
extensionsFull["topJScore"] = extensionsFull.allJHitsWithScore.str.extract("^[^(]+\((?P<extended>[0-9.]+)\)", expand=False).apply(pd.to_numeric)

extensionsLEx = extensionsFull[extensionsFull.targetDescriptions.str.contains("LExtended")]
extensionsLEx = pd.concat([extensionsLEx, extensionsLEx.targetDescriptions.str.extract("LExtended\((?P<extended>[0-9]+)\)", expand=True).apply(pd.to_numeric)], axis=1)
extensionsLEx["actualExtension"] = extensionsLEx.apply(lambda row: row.readSequence[:row.extended], axis=1)
extensionsLEx["expectedExtension"] = extensionsLEx.apply(lambda row: row.trueCDR3[:row.extended], axis=1)
# extensionsLEx[extensionsLEx.actualExtension != extensionsLEx.expectedExtension]

extensionsREx = extensionsFull[extensionsFull.targetDescriptions.str.contains("[RM]Extended")]
extensionsREx = pd.concat([extensionsREx, extensionsREx.targetDescriptions.str.extract("\[(?P<offset>[0-9]+)\] \+ [RM]Extended\((?P<extended>[0-9]+)\)", expand=True).apply(pd.to_numeric)], axis=1)
extensionsREx["actualExtension"] = extensionsREx.apply(lambda row: row.readSequence[row.offset:(row.offset + row.extended)], axis=1)
extensionsREx["offsetInTrueCDR3"] = extensionsREx.offset - extensionsREx.CDR3End + extensionsREx.trueCDR3.str.len()
extensionsREx["offsetInTrueCDR3L"] = extensionsREx.offset - extensionsREx.CDR3Begin
extensionsREx.loc[extensionsREx.offsetInTrueCDR3.isnull(), "offsetInTrueCDR3"] = extensionsREx.loc[extensionsREx.offsetInTrueCDR3.isnull(), "offsetInTrueCDR3L"] 
extensionsREx["expectedExtension"] = extensionsREx.apply(lambda row: row.trueCDR3[int(row.offsetInTrueCDR3):int(row.offsetInTrueCDR3 + row.extended)], axis=1)
# extensionsREx[extensionsREx.actualExtension != extensionsREx.expectedExtension]

clones = 0
for line in open(assembleReport,'r'):
    m = re.search("Reads used in clonotypes, percent of total: ([0-9]+)", line)
    if m:
        clones = int(m.group(1))
        break

result = {
    "inputFile": inputFile,
    "clonesTotal": clones,
    "totalRExtensions": len(extensionsREx),
    "totalLExtensions": len(extensionsLEx),
    "falseRExtensions": int((extensionsREx.actualExtension != extensionsREx.expectedExtension).sum()),
    "falseLExtensions": int((extensionsLEx.actualExtension != extensionsLEx.expectedExtension).sum())
}

print(json.dumps(result))
