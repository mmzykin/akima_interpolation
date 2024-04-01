CUNITS = akima.c
all: $(CUNITS)
	gcc $(CUNITS) -lm  -o akima
clean:
	rm akima *_out

