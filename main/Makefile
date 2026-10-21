SRC      = $(wildcard *.cpp)
OBJ      = $(SRC:.cpp=.o)
DEP      = $(patsubst %.cpp,%.d,$(SRC))
NAME     = $(SRC:.cpp=)

ACCURT_KERNEL_DIR ?= ../external/accurt-kernel
ACCURT_KERNEL_AVAILABLE := $(if $(wildcard $(ACCURT_KERNEL_DIR)/kernel.f95),1,)
KERNEL_CONFIG := .accurt-kernel-$(if $(ACCURT_KERNEL_AVAILABLE),enabled,disabled)

ifneq ($(ACCURT_KERNEL_AVAILABLE),)
FC                 := gfortran
ACCURT_KERNEL_OBJ  := ../accurt_kernel_api/kernel.o
FORTRAN_LIB_DIR    := $(dir $(shell $(FC) -print-file-name=libgfortran.dylib))
LDLIBS             += -L$(FORTRAN_LIB_DIR) -lgfortran -lquadmath -framework Accelerate
CPPFLAGS           += -DFLICK_HAS_ACCURT_KERNEL
endif

all:	obj link test

clean:
	@rm -f $(OBJ) $(DEP) ./$(NAME) ./*~ ./obj ./#*# .accurt-kernel-enabled .accurt-kernel-disabled
link: $(ACCURT_KERNEL_OBJ)
	$(FLICK_COMPILER) -o $(NAME) $(OBJ) $(ACCURT_KERNEL_OBJ) $(LDLIBS)

obj:	$(OBJ) $(ACCURT_KERNEL_OBJ)

$(OBJ): $(KERNEL_CONFIG)

$(KERNEL_CONFIG):
	@rm -f .accurt-kernel-enabled .accurt-kernel-disabled
	@touch $@

$(ACCURT_KERNEL_OBJ):
	$(MAKE) -C ../accurt_kernel_api kernel.o

test:
	@./$(NAME)

-include $(DEP) 

%.o: %.cpp Makefile
	$(FLICK_COMPILER) $(CPPFLAGS) -MMD -MP -c $< -o $@

