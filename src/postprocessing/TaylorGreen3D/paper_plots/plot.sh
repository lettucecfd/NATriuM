#! /bin/sh
files="15000_dissipation 1600_dissipation 800_dissipation 200_dissipation dissipation eff_Re"
for filename in $files
do
  R --no-save <  ${filename}.R
  pdflatex ${filename}.tex
  evince ${filename}.pdf &
done

