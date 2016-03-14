cppname = Project1
# compilecpp = g++
compilecpp = gfortran
# flags=-Wall -Wextra -Wconversion -std=c++11 -O3 -o
flags=-cpp -Wall -Wextra -Wconversion -O3 -fdefault-real-8 -o

all: compilecpp run

compilecpp:
	$(compilecpp) $(cppname).f90 $(flags) $(cppname).exec
# 	$(compilecpp) $(cppname).cpp $(flags) $(cppname).exec

run:
	./$(cppname).exec

plot:
	python *.py
	gnome-open *.pdf

clean:
	rm *.exec 
