


COMPILER=gcc
OPTIONS=-pedantic -Wall -Wextra -Werror -Wshadow -Wconversion -Wunreachable-code
COMPILE=$(COMPILER) $(OPTIONS)


kalman_filter: kalman_filter.cpp sensor_handlers.cpp math_util.cpp
	$(COMPILE) $^ -o $@ -lm


clean:
	rm kalman_filter

.PHONY: clean

