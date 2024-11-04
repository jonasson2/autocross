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

# Detect the OS
ifeq ($(OS),Windows_NT)
    RM = del /F /Q
    RMDIR = rmdir /S /Q
    SEP = /
    EXE = .exe
else
    RM = rm -f
    RMDIR = rm -rf
    SEP = /
    EXE =
endif

SOURCES = $(wildcard $(SRC)$(SEP)*.f90)
$(info $(SEP))
$(info $(SOURCES))
$(info $(SRC))
OBJECTS = $(patsubst $(SRC)$(SEP)%.f90,$(OBJ)$(SEP)%.o,$(SOURCES))
EXECUTABLES = $(BIN)$(SEP)pearsont3$(EXE) $(BIN)$(SEP)redfit-x$(EXE)
$(info $(OBJECTS))

all: $(EXECUTABLES)

$(OBJ)$(SEP)%.o: $(SRC)$(SEP)%.f90
	-mkdir $(OBJ)
	-mkdir $(MOD)
	$(FC) $(CFLAGS) $(FLAGS) -c $< -o $@

$(BIN)$(SEP)pearsont3$(EXE): $(OBJECTS)
	-mkdir $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

$(BIN)$(SEP)redfit-x$(EXE): $(OBJECTS)
	-mkdir $(BIN)
	$(FC) $(LFLAGS) $(FLAGS) $^ -o $@

pears:
	cd run && time pearsont3 && cd ..

clean:
	-$(RM) $(OBJ)$(SEP)*.o $(MOD)$(SEP)*.mod $(BIN)$(SEP)*$(EXE)
	-$(RMDIR) $(OBJ) $(MOD)

.PHONY: all clean

