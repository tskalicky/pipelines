#!/bin/bash
## Workflow in this folder
## Annotate all samples and intersect foreground with
## selection (bayesian statistics - done before).
#########################################################

## 1. Annotate all samples
##    - the annotation goes:
##      tRNA -> miRNA -> snRNA -> rRNA -> miscRNA -> snoRNA -> mRNA -> lincRNA
##    - tRNA has offset 30nt!
##    - mRNA = CDS + UTR

## Annotate
echo "Annotating files"
for I in mut1 mut2 mut3; do
  ./annotate.sh $I foreground.bed.gz &
done
wait

for I in mut1 mut2 mut3; do
  ./annotate.sh $I foreground.uadd.bed.gz &
done
wait

## Intersect reads with selected significant regions
echo "Intersecting with selected regions"
for I in mut1 mut2 mut3; do
  ./select_foreground.sh $I &
done
wait

## Correct the multiplicity
## POZN. BedTools se chová divně = pokud intersect vícekrát, nahlásí výsledek vícekrát.
## proto je potřeba normalizace na anotaci
echo "Counting multiplicity"
for I in mut1 mut2 mut3; do
  ./count_multiplicity.pl $I foreground.ann.bed.gz &
  ./count_multiplicity.pl $I foreground.ann.sel.bed.gz &
  ./count_multiplicity.pl $I foreground.uadd.ann.bed.gz &
  ./count_multiplicity.pl $I foreground.uadd.ann.sel.bed.gz &
  wait
done

## Count multiplicity for transcripts
for I in mut1 mut2 mut3; do
  ./count_multiplicity.pl $I foreground.ann_transcripts.bed.gz &
  ./count_multiplicity.pl $I foreground.uadd.ann_transcripts.bed.gz &
  wait
done
