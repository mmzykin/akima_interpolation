CUNITS = akima.c
all: $(CUNITS)
	gcc $(CUNITS) -lm  -o akima
clean:
	rm -f akima *_out
test: clean debug
	./akima -f input -s 1 
debug: $(CUNITS)
	gcc $(CUNITS) -lm -D_DEBUG -g -o akima

