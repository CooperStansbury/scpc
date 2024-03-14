#!/bin/bash

# a script to annotate aligned monomers
# based on https://github.com/epi2me-labs/pore-c-py

BAM=""
THREADS=4
OUTPUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
     -b|--bam)
      shift
      BAM="$1"
      ;; 
    -t|--threads)
      shift
      THREADS="$1"
      ;;
    -o|--output)
      shift
      OUTPUT="$1"
      ;;
    *)
      echo "Invalid option: $1"
      exit 1
      ;;
  esac
  shift
done

echo "Input BAM:"  "${BAM}"
echo "Threads:" "${THREADS}"
echo "Output:" "${OUTPUT}"


pore-c-py annotate "${BAM}" "${OUTPUT}" --force --monomers --stdout true --summary --chromunity > "${OUTPUT}".ns.bam