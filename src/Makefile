


COMPILER=gcc
OPTIONS=-pedantic -Wall -Wextra -Werror -Wshadow -Wconversion -Wunreachable-code
COMPILE=$(COMPILER) $(OPTIONS)


all: kalman_filter testmath

kalman_filter: kalman_filter.c sensor_handlers.c math_util.c
	$(COMPILE) $^ -o $@ -lm

testmath: testmath.c math_util.c
	$(COMPILE) $^ -o $@ -lm

clean:
	rm kalman_filter testmath

.PHONY: clean

