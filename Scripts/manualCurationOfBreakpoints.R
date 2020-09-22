setwd("/")
args = commandArgs(trailingOnly=TRUE)
#args= c("/Users/zeidh/Desktop/BreakpointerGithub/Output/Breakpoints/","Aug28_breakpoints.txt")
bind <- read.table(paste0(args[1],args[2]),header=T)

#ALL
#countedSCEs <- bind[-c(2:3,5:7,12,17:18,29,34:35,37,42,46,47,50,66,71,74,87,109,124,128,129,138,149:151,209,240,247,350,383,418,419,464:467,484,495,496,517,522,539,546,547,581,595:597,629,632,634,636:637,650,651,658,659,673,675,676,691,692,708,709,756,757,770,771,771,789,801,806,807,827,828,836,837,835,839,840,852,853,857,861,862,864,865,867,868,877,878,882:885,887,888:891,912,913,917:927,941:944,959:966,975,977:980,983,996:1001,1007,1020,1021,1058,1080,1085,1086,1093,1094,1101,1102,1122:1125,1159,1171,1172,1176,1192,1193,1179,1198:1201,1205,1209,1215,1217,1218,1222:1224),]

#Aug28
countedSCEs <- bind[-c(3,22:24,27,33,35,118,123:126,178,182,191,224,226,251,279,287,293,302:308,320,321,334,342,348,362,395,399,400,402,405,412,443,448),]
rows <- c()
#removing chr15 and chr16 inversions at 60-90Mb and 21-23Mb
#removes 258 inversion calls (129) from 1062 leaving behind 804 SCEs
for (row in 1:nrow(countedSCEs)){
	if (countedSCEs[row,1]=="chr15"){
		if (countedSCEs[row,2] <= 60000000 & countedSCEs[row,3] >= 90000000 ){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] >= 60000000 & countedSCEs[row,2] <= 90000000 & countedSCEs[row,3] <= 90000000 & countedSCEs[row,3] >= 60000000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] >= 60000000 & countedSCEs[row,2] <= 90000000 & countedSCEs[row,3] >= 90000000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] <= 60000000 & countedSCEs[row,3] <= 90000000 & countedSCEs[row,3] >= 60000000){
			rows <- append(row,rows)
		}
	}
	else if (countedSCEs[row,1]=="chr16"){
		if (countedSCEs[row,2] <= 21000000 & countedSCEs[row,3] >= 23000000 ){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] >= 21000000 & countedSCEs[row,2] <= 23000000 & countedSCEs[row,3] <= 23000000 & countedSCEs[row,3] >= 21000000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] >= 21000000 & countedSCEs[row,2] <= 23000000 & countedSCEs[row,3] >= 23000000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] <= 21000000 & countedSCEs[row,3] <= 23000000 & countedSCEs[row,3] >= 21000000){
			rows <- append(row,rows)
		}
	}
	else if (countedSCEs[row,1]=="chr8"){
		if (countedSCEs[row,2] >= 70569000 & countedSCEs[row,2] <= 70800000 & countedSCEs[row,3] <= 70800000 & countedSCEs[row,3] >= 70569000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] >= 70569000 & countedSCEs[row,2] <= 70800000 & countedSCEs[row,3] >= 70800000){
			rows <- append(row,rows)
		}
		else if (countedSCEs[row,2] <= 70569000 & countedSCEs[row,3] <= 70800000 & countedSCEs[row,3] >= 70569000){
			rows <- append(row,rows)
		}
	}
}
countedSCEs <- (countedSCEs[-c(rows),])


