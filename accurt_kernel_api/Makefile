CPP_SRC  := $(wildcard *.cpp)
CPP_OBJ  := $(CPP_SRC:.cpp=.o)
DEP      := $(CPP_SRC:.cpp=.d)
NAME     := $(CPP_SRC:.cpp=)

CXX      := $(if $(FLICK_COMPILER),$(FLICK_COMPILER),clang++)
FC       := gfortran
CXXFLAGS ?= -std=c++20
FFLAGS   += -O2 -finit-integer=0

ACCURT_KERNEL_DIR ?= ../external/accurt-kernel
FORTRAN_SRC       := $(ACCURT_KERNEL_DIR)/kernel.f95
ACCURT_KERNEL_AVAILABLE := $(if $(wildcard $(FORTRAN_SRC)),1,)
KERNEL_CONFIG     := .accurt-kernel-$(if $(ACCURT_KERNEL_AVAILABLE),enabled,disabled)

ifneq ($(ACCURT_KERNEL_AVAILABLE),)
FORTRAN_OBJ       := kernel.o
ifeq ($(shell uname -s),Linux)
LDFLAGS           += -O2
LDLIBS            += -lgfortran -lquadmath -llapack -lblas
else
FORTRAN_LIB_DIR   := $(dir $(shell $(FC) -print-file-name=libgfortran.dylib))
LDLIBS            += -L$(FORTRAN_LIB_DIR) -lgfortran -lquadmath -framework Accelerate
endif
CXXFLAGS          += -DFLICK_HAS_ACCURT_KERNEL
endif

.PHONY: all clean obj link test

all: test

clean:
	@rm -f $(CPP_OBJ) $(FORTRAN_OBJ) $(DEP) $(NAME) *.mod *~ obj ./#*# tmp.txt .accurt-kernel-enabled .accurt-kernel-disabled
	@rm -fr tmp tmpOutput/ tmpMaterials/ flick_tmp/
	@rm -f toa_reflectance_configuration

link: $(CPP_OBJ) $(FORTRAN_OBJ)
	$(CXX) $(LDFLAGS) -o $(NAME) $(CPP_OBJ) $(FORTRAN_OBJ) $(LDLIBS)

obj:	$(CPP_OBJ) $(FORTRAN_OBJ)

$(CPP_OBJ): $(KERNEL_CONFIG)

$(KERNEL_CONFIG):
	@rm -f .accurt-kernel-enabled .accurt-kernel-disabled
	@touch $@

test: link
	@if [ -z "$(ACCURT_KERNEL_AVAILABLE)" ]; then \
		echo "Skipping accurt_kernel_api tests: $(FORTRAN_SRC) not found."; \
	else \
		./$(NAME); \
	fi

-include $(DEP) 

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(FORTRAN_OBJ): $(FORTRAN_SRC) Makefile
	$(FC) $(FFLAGS) -c $< -o $@

