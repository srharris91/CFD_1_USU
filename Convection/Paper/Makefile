# type "make" command in Unix to create asme2e.pdf file 
all: file clean

file:
	pdflatex Project2
	bibtex Project2
	pdflatex Project2
	pdflatex Project2
	# dvips -o Project2.ps Project2
	# ps2pdf Project2.ps Project2.pdf
	gnome-open Project2.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
