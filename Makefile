CC = gcc
CXXFLAGS = -Wall -O3
LDFLAGS =
LIBS = -lm -lgsl -lgslcblas

SRC_FILES = vector.c kurtic.c solver.c
OBJS = $(SRC_FILES:%.c=%.o)

NAME = excess

%.o: %.c %.h
	$(CC) $(CXXFLAGS) -c $< -o $@

$(NAME): main.c $(OBJS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

tests: test.c $(OBJS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

excess_calculator: excess_calculator.c $(OBJS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRC_FILES)
	$(CC) -MM $^ > $@

clean:
	rm -f *.o
	rm -f deps.mk
	rm -f $(NAME)
	rm -f tests
