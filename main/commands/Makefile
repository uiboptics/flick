SRC      = $(wildcard *.cpp)
OBJ      = $(SRC:.cpp=.o)
DEP      = $(patsubst %.cpp,%.d,$(SRC))
NAME     = $(SRC:.cpp=)

.PHONY: all clean parent obj link test

all:	test

clean:
	@rm -f $(OBJ) $(DEP) ./$(NAME) ./*~ ./obj
parent:
	$(MAKE) -C .. obj link

link:	obj
	$(FLICK_COMPILER) -o $(NAME) $(OBJ)

obj:	$(OBJ)

test:	parent link
	@./$(NAME)

-include $(DEP) 

%.o: %.cpp Makefile
	$(FLICK_COMPILER) -MMD -MP -c $< -o $@

