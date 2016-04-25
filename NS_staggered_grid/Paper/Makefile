# type "make" command in Unix to create asme2e.pdf file 
all: file clean

file:
	pdflatex Project3
	bibtex Project3
	pdflatex Project3
	pdflatex Project3
	# dvips -o Project3.ps Project3
	# ps2pdf Project3.ps Project3.pdf
	gnome-open Project3.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
