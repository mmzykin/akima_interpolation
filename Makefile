CUNITS = akima.c
all: $(CUNITS)
	gcc $(CUNITS) -lm -Wall -Wextra -Wpedantic -fsanitize=address -g -o akima
clean:
	rm akima

