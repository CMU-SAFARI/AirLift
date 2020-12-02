CC=gcc
CFLAGS = -c -O3 -Wall -msse -msse2 
LDFLAGS = -lz -lm 
SOURCES = baseFAST.c CommandLineParser.c Common.c HashTable.c MrFAST.c Output.c Reads.c RefGenome.c GreedySnake.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = mrfast


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* mrfast
