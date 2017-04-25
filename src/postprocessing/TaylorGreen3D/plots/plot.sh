python calculate_dEdt_inc.py
R --no-save < plot_kE_keq.R
pdflatex dEdt_keq.tex
evince dEdt_keq.pdf
convert -density 400 dEdt_keq.pdf dEdt_keq.png
