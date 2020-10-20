#!/usr/bin/awk -f

BEGIN {
    OFS=","
	print "time", "id", "type", "x", "y", "z", "ix", "iy", "iz"
}

{
if        ($0 == "ITEM: TIMESTEP") {
	getline timestep
	inCoords = 0
} else if (inCoords) {
	print timestep, $1, $2, $3, $4, $5, $6, $7, $8
} else if ($0 ~ "ITEM: ATOMS") {
	inCoords = 1
}
}
