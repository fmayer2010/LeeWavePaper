MAIN = nondimensional

all:	dvi

dvi:
	latex $(MAIN)

bib:
	latex $(MAIN)
	bibtex $(MAIN)
	latex $(MAIN)
	latex $(MAIN)

ps:     dvi
	dvips -Ppdf -G0 $(MAIN) -o

pdf:    ps
	ps2pdf $(MAIN).ps $(MAIN).pdf

view:	dvi
	xdvi $(MAIN) &

clean:
	rm -f *.aux *.bbl *.blg *.log *~ \#*\# *.toc *.dvi *.ps 

clobber:	clean
	rm -f *.pdf
