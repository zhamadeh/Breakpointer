#!/bin/bash

#take metadata file with libraries characterized as g (good),p (poor), d (shallow) and move to subdir
Rscript Scripts/createSubDirOfGoodLibraries.R Input/metrics/Aug28.txt  Input/RData_unfiltered/ Input/RData_good/

Rscript Scripts/blacklistCentromeresFromRData.R Input/RData_good/ Input/RData_blacklisted/ Output/BPR_breaksPlots/

Rscript Scripts/collectBreakpoints.R Input/RData_blacklisted/ Output/Breakpoints/ Input/Metrics/ Aug28_breakpoints.txt

Rscript Scripts/plotting.R Output/Breakpoints/Aug28_breakpoints.txt

Rscript Scripts/savingBreakpointsToEnrichmentDir.R Output/Breakpoints/Aug28_breakpoints.txt
