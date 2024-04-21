FC = gfortran
SRC = src
OBJ = obj
BIN = bin
MOD = modules

CFLAGS = -O3 -J$(MOD)
LFLAGS = -I$(MOD)

SOURCES = $(wildcard $(SRC)/*.f90)
OBJECTS = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SOURCES))
EXECUTABLES = $(BIN)/pearsont3 $(BIN)/redfit-x

all: $(EXECUTABLES)

$(BIN)/pearsont3: $(OBJ)/modules.o $(OBJ)/pearsont3.o
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $^ -o $@

$(BIN)/redfit-x: $(OBJ)/modules.o $(OBJ)/redfit-x.o
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.f90
	mkdir -p $(OBJ) $(MOD)
	$(FC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ)/* $(MOD)/* $(BIN)/*

.PHONY: all clean
