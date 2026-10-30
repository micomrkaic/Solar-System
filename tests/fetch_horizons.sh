#!/bin/sh
# Fetch heliocentric ecliptic-J2000 position vectors for the 8 planets from
# JPL Horizons at one epoch and write tests/horizons.csv for `make test`.
#
#   tests/fetch_horizons.sh [JD_TDB]      default 2458849.5 = 2020-01-01 00:00
#
# The reference frame matches the simulation: Sun-centred (500@10),
# ecliptic and mean equinox of J2000, AU.  TDB vs TT differs by < 2 ms.
# Bodies are the planetary-system barycentres 1–8 (3 = Earth–Moon
# barycentre), matching the system masses the simulation integrates.
# With horizons_j2000.h present (tests/fetch_j2000_state.sh) the check is
# strict; from the approximate mean-elements table it is informational.
set -eu
JD="${1:-2458849.5}"
OUT="$(dirname "$0")/horizons.csv"
API="https://ssd.jpl.nasa.gov/api/horizons.api"
echo "# name,JD_TDB,X_AU,Y_AU,Z_AU   (JPL Horizons, heliocentric ecliptic J2000)" > "$OUT"
i=0
for name in Mercury Venus Earth Mars Jupiter Saturn Uranus Neptune; do
    i=$((i+1))
    id="$i"
    line=$(curl -sS "$API?format=text&COMMAND='$id'&EPHEM_TYPE='VECTORS'&CENTER='500@10'&TLIST='$JD'&REF_PLANE='ECLIPTIC'&VEC_TABLE='1'&OUT_UNITS='AU-D'&CSV_FORMAT='YES'" \
           | awk '/\$\$SOE/{f=1;next} /\$\$EOE/{f=0} f' | head -1)
    [ -n "$line" ] || { echo "no data for $name ($id)" >&2; exit 1; }
    # CSV: JDTDB, Calendar Date (TDB), X, Y, Z,
    echo "$line" | awk -F',' -v n="$name" '{gsub(/ /,"",$1);gsub(/ /,"",$3);gsub(/ /,"",$4);gsub(/ /,"",$5); printf "%s,%s,%s,%s,%s\n", n, $1, $3, $4, $5}' >> "$OUT"
done
echo "wrote $OUT:"; cat "$OUT"
