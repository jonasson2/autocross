# Breytt af KBOx
FC = gfortran
SRC = src
OBJ = obj
BIN = .
MOD = modules
MODOBJ = $(addprefix $(OBJ)/,pearsont3mod.o modules.o pearsont3sub_mod.o)

CFLAGS = -J$(MOD)
LFLAGS = -I$(MOD)

# Add debug flags for debugging
DEBUGFLAGS = -g -O0 -fbacktrace -fcheck=all
RELEASEFLAGS = -O3
FLAGS = $(if $(DEBUG),$(DEBUGFLAGS),$(RELEASEFLAGS))

SOURCES = $(wildcard $(SRC)/*.f90)
OBJECTS = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(SOURCES))
EXECUTABLES = $(BIN)/pearsont3 $(BIN)/redfit-x
DYNAMICLIBRARIES=pearsont3.so
# $(info $(OBJECTS))

all: $(EXECUTABLES) $(DYNAMICLIBRARIES)

$(OBJ)/pearsont3.o: $(MODOBJ)

pearsont3.so: $(OBJECTS)
	$(FC) $(LFLAGS) -shared -fPIC $(FLAGS) $^ -o $@  

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
	rm -rf $(OBJ)/* $(MOD)/* $(EXECUTABLES) $(DYNAMICLIBRARIES)

pearson:
	cd run && time pearsont3 && cd ..

.PHONY: all clean
