.PHONY: clean

GCC = g++ 
CFLAGS = -O3 -mcmodel=medium -ffast-math -Wall
LDFLAGS = -I.
objs = main.o


rocket: $(objs)
	$(GCC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp 
	$(GCC) $(CFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -f *.o *~
