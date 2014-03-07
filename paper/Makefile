pdflatex: *.tex
	rm -f document.aux document.pdf document.log document.bbl
	pdflatex document > pdflatex.err
	bibtex document.aux
	pdflatex document >>pdflatex.err
	pdflatex document >>pdflatex.err
	pdflatex document >>pdflatex.err

all:
	pdflatex

clean:
	rm -f document.aux document.pdf document.log document.bbl

supp: supplement.tex *.eps
	rm -f supplement.pdf supplement.aux *-converted-to.pdf
	pdflatex supplement.tex
	bibtex supplement.aux
	pdflatex supplement
	pdflatex supplement
