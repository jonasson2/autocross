# Breytt af KBOx
FC = gfortran
SRC = src
OBJ = obj
BIN = .
MOD = modules
MODOBJ = $(addprefix $(OBJ)/,pearsont3mod.o modules.o pearsont3sub_mod.o)

CFLAGS = -O3 -J$(MOD)
LFLAGS = -I$(MOD)

# Add debug flags for debugging
DEBUGFLAGS = -g -fbacktrace -fcheck=all
RELEASEFLAGS = -O3
FLAGS = $(if $(DEBUG),$(DEBUGFLAGS),$(RELEASEFLAGS))

SOURCES = $(wildcard $(SRC)/*.f90)
OBJECTS = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SOURCES))
EXECUTABLES = $(BIN)/pearsont3 $(BIN)/redfit-x
$(info $(OBJECTS))

all: $(EXECUTABLES)

$(OBJ)/pearsont3.o: $(MODOBJ)

$(BIN)/pearsont3: $(OBJECTS) 
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(BIN)/redfit-x: $(OBJECTS)
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(OBJ)/%.o: $(SRC)/%.f90
	mkdir -p $(OBJ) $(MOD)
	$(FC) $(CFLAGS) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ)/* $(MOD)/* $(BIN)/*.exe

pears:
	cd run && time pearsont3 && cd ..

.PHONY: all clean

