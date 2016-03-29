code = Project3
#module = Project2_mod.f90
compiler = gfortran
# flags=-Wall -Wextra -Wconversion -std=c++11 -O3 -o
# flags=-cpp -Wall -Wextra -Wconversion -J./output/ -fopenmp -O3 -fdefault-real-8 -o
flags=-cpp -Wall -Wextra -Wconversion -J./output/ -O3 -fdefault-real-8
output = output/
EXEC = Project3.exec
OBJS = Types.o Project3.o

all: exec run plot

compile:
	$(compiler) $(module) $(code).f90 $(flags) -o $(output)$(code).exec
# 	$(compiler) $(code).cpp $(flags) $(code).exec
exec: $(OBJS)
	@echo ''
	@echo 'Linking $(output)$(EXEC) from $(OBJS)'
	@echo ''
	@$(compiler) $(flags) $(OBJS) -o $(output)$(EXEC)

run:
	./$(output)$(code).exec

plot:
	python *.py

clean:
	rm -f *.exec *.txt *.pdf *.mod output/* *.o




.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	@echo 'Compiling [$@] from [$<] using $(comp) $(flags)'
	@$(compiler) $(flags) $*.f90 -c -o $*.o
