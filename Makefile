# Breytt af KBOx
FC = gfortran
SRC = src
OBJ = obj
BIN = bin
MOD = modules

CFLAGS = -O3 -J$(MOD)
LFLAGS = -I$(MOD)
FLAGS =

SOURCES = $(wildcard $(SRC)/*.f90)
OBJECTS = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SOURCES))
EXECUTABLES = $(BIN)/pearsont3 $(BIN)/redfit-x

all: $(EXECUTABLES)

$(BIN)/pearsont3: $(OBJ)/modules.o $(OBJ)/pearsont3.o
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(BIN)/redfit-x: $(OBJ)/modules.o $(OBJ)/redfit-x.o
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.f90
	mkdir -p $(OBJ) $(MOD)
	$(FC) $(CFLAGS) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ)/* $(MOD)/* $(BIN)/*

pears:
	cd run && time pearsont3 && cd ..

.PHONY: all clean
