# type "make" command in Unix to create asme2e.pdf file 
all: file clean

file:
	pdflatex Project1
	bibtex Project1
	pdflatex Project1
	pdflatex Project1
	# dvips -o Project1.ps Project1
	# ps2pdf Project1.ps Project1.pdf
	gnome-open Project1.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
