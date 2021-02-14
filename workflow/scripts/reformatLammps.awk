#!/usr/bin/awk -f

BEGIN {
    OFS=","
    N = 0
    outRep = 0
    if (nSteps == "") {
        nSteps = 10
    }
    if (prefix == "") {
        print "Must set variable: prefix" > "/dev/stderr"
        exit 1
    }
}

{
if        ($0 == "ITEM: TIMESTEP") {
	getline timestep
	inCoords = 0
    if (N % nSteps == 0) {
        out = prefix""outRep".csv.gz"
        print "timestep", "ID", "type", "x", "y", "z", "ix", "iy", "iz" | "gzip >" out
        outRep += 1
    }
    N += 1
} else if (inCoords) {
	print timestep, $1, $2, $3, $4, $5, $6, $7, $8 | "gzip >" out
} else if ($0 ~ "ITEM: ATOMS") {
	inCoords = 1
}
}
