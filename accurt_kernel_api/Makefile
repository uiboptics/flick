CPP_SRC  := $(wildcard *.cpp)
CPP_OBJ  := $(CPP_SRC:.cpp=.o)
DEP      := $(CPP_SRC:.cpp=.d)
NAME     := $(CPP_SRC:.cpp=)

CXX      := $(if $(FLICK_COMPILER),$(FLICK_COMPILER),clang++)
FC       := gfortran
CXXFLAGS ?= -std=c++20
FFLAGS   += -O3 -fallow-argument-mismatch -finit-integer=0

ACCURT_KERNEL_DIR ?= /Users/borgehamre/f/external/accurt-kernel
FORTRAN_SRC       := $(ACCURT_KERNEL_DIR)/kernel.f95
FORTRAN_OBJ       := kernel.o
FORTRAN_LIB_DIR   := $(dir $(shell $(FC) -print-file-name=libgfortran.dylib))
LDLIBS            += -L$(FORTRAN_LIB_DIR) -lgfortran -lquadmath -framework Accelerate

.PHONY: all clean obj link test

all: test

clean:
	@rm -f $(CPP_OBJ) $(FORTRAN_OBJ) $(DEP) $(NAME) *.mod *~ obj ./#*# tmp.txt
	@rm -fr tmp tmpOutput/ tmpMaterials/ flick_tmp/
	@rm -f toa_reflectance_configuration

link: $(CPP_OBJ) $(FORTRAN_OBJ)
	$(CXX) $(LDFLAGS) -o $(NAME) $(CPP_OBJ) $(FORTRAN_OBJ) $(LDLIBS)

obj:	$(CPP_OBJ) $(FORTRAN_OBJ)

test: link
	@./$(NAME)

-include $(DEP) 

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(FORTRAN_OBJ): $(FORTRAN_SRC) Makefile
	$(FC) $(FFLAGS) -c $< -o $@

