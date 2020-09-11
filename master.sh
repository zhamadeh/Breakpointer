#!/bin/bash

#take metadata file with libraries characterized as g (good),p (poor), d (shallow) and move to subdir
Rscript Scripts/createSubDirOfGoodLibraries.R /Users/zeidh/Desktop/Data/METRICS/Aug18.txt  /Users/zeidh/Desktop/Data/RDATA/Aug18-2020/unfiltered_fromServer/ /Users/zeidh/Desktop/Data/RDATA/Aug18-2020/good/

cp /Users/zeidh/Desktop/Data/RDATA/Aug18-2020/good/* /Users/zeidh/Desktop/Data/RDATA/ALL/

Rscript Scripts/blacklistCentromeresFromRData.R /Users/zeidh/Desktop/Data/RDATA/ALL/ /Users/zeidh/Desktop/Data/RDATA/ALL_BL/ /Users/zeidh/Desktop/Data/PLOTS/

Rscript Scripts/collectBreakpoints.R /Users/zeidh/Desktop/Data/RDATA/ALL_BL /Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/ /Users/zeidh/Desktop/Data/METRICS
