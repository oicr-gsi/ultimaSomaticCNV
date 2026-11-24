#!/bin/bash
# Summarise a provisioned output directory into a form a baseline can be compared against:
# the file inventory, the VCF record count per FILTER value and each BED's line count.
# No sizes, paths, headers, image contents or floating-point values, so the baseline does
# not drift.
set -o nounset
set -o errexit
set -o pipefail
export LC_ALL=C

cd "$1" || exit 1

find . -type f | sed 's|^\./||' | sort

for vcf in *.vcf.gz; do
    [ -e "$vcf" ] || continue
    echo "## ${vcf}"
    zcat "$vcf" | grep -v '^#' | awk -F'\t' '{ print $7 }' | sort | uniq -c | sort -k2,2
done

for bed in *.bed; do
    [ -e "$bed" ] || continue
    echo "## ${bed} $(wc -l < "$bed") lines"
done
