#!/bin/bash
## bash select_foreground.sh mut1 notu

homeSmp="$HOME/Projects/dis3l2/samples"
homeSel="$HOME/Projects/dis3l2/statistics/mihaela/selection/output"
smp=$1   # Sample to analyze: mut1

if [ $# -ne 1 ]; then
  echo "Please set one argument: sample_to_analyze"
  exit 0
fi

## Intersect regular annotation
intersectBed -wa -s \
     -a "$homeSmp/$smp/foreground.uadd.ann.bed.gz" \
     -b "$homeSel/selected.bed.gz" \
      | gzip \
      > "$homeSmp/$smp/foreground.uadd.ann.sel.bed.gz"

intersectBed -wa -s \
     -a "$homeSmp/$smp/foreground.ann.bed.gz" \
     -b "$homeSel/selected.bed.gz" \
      | gzip \
      > "$homeSmp/$smp/foreground.ann.sel.bed.gz"

## Intersect for transcript annotation
intersectBed -wa -s \
     -a "$homeSmp/$smp/foreground.uadd.ann_transcripts.bed.gz" \
     -b "$homeSel/selected.bed.gz" \
      | gzip \
      > "$homeSmp/$smp/foreground.uadd.ann_transcripts.sel.bed.gz"

intersectBed -wa -s \
     -a "$homeSmp/$smp/foreground.ann_transcripts.bed.gz" \
     -b "$homeSel/selected.bed.gz" \
      | gzip \
      > "$homeSmp/$smp/foreground.uadd.ann_transcripts.sel.bed.gz"
