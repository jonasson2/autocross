# Breytt af KBOx
FC = gfortran
SRC = src
OBJ = obj
BIN = .
MOD = modules

CFLAGS = -J$(MOD)
LFLAGS = -I$(MOD)

# Add debug flags for debugging
DEBUGFLAGS = -g -O0 -fbacktrace -fcheck=all
RELEASEFLAGS = -O3
FLAGS = $(if $(DEBUG),$(DEBUGFLAGS),$(RELEASEFLAGS))

P3_SOURCES = pearsont3.f90 p3_subroutine.f90 p3_modules.f90 common_modules.f90
RX_SOURCES = redfitx.f90 rx_modules.f90 common_modules.f90

P3_OBJECTS = $(patsubst %.f90,$(OBJ)/%.o,$(P3_SOURCES))
RX_OBJECTS = $(patsubst %.f90,$(OBJ)/%.o,$(RX_SOURCES))
OBJECTS = $(sort $(P3_OBJECTS) $(RX_OBJECTS))

EXECUTABLES = $(BIN)/pearsont3 $(BIN)/redfit-x
DYNAMICLIBRARIES=pearsont3.so

all: $(EXECUTABLES) $(DYNAMICLIBRARIES)

$(OBJ)/%.o: $(SRC)/%.f90
	mkdir -p $(OBJ) $(MOD)
	$(FC) $(CFLAGS) $(FLAGS) -MMD -MP -c $< -o $@

pearsont3.so: $(filter-out $(OBJ)/pearsont3.o, $(P3_OBJECTS))
	$(FC) $(LFLAGS) -shared -fPIC $(FLAGS) $^ -o $@  

$(BIN)/pearsont3: $(P3_OBJECTS)
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(BIN)/redfit-x: $(RX_OBJECTS)
	mkdir -p $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

clean:
	rm -rf $(OBJ)/* $(MOD)/* $(EXECUTABLES) $(DYNAMICLIBRARIES)

pearson:
	cd run && time pearsont3 && cd ..

.PHONY: all clean
