code = Project3
#module = Project2_mod.f90
compiler = gfortran
# flags=-Wall -Wextra -Wconversion -std=c++11 -O3 -o
# flags=-cpp -Wall -Wextra -Wconversion -J./output/ -fopenmp -O3 -fdefault-real-8 -o
flags=-cpp -Wall -Wextra -Wconversion -J./output/ -O3 -fdefault-real-8 -o
output = output/

all: compile run

compile:
	$(compiler) $(module) $(code).f90 $(flags) $(output)$(code).exec
# 	$(compiler) $(code).cpp $(flags) $(code).exec

run:
	./$(output)$(code).exec

plot:
	python *.py

clean:
	rm -f *.exec *.txt *.pdf *.mod output/*
