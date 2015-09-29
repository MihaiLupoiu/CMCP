CC = g++
LFLAGS = -Wall -lpthread
CFLAGS = -Wall -c

ejer2 : ejer2.o
	$(CC) $(LFLAGS) ejer2.o -o ejer2

ejer3 : ejer3.o
	$(CC) $(LFLAGS) ejer3.o -o ejer3

ejer2.o :
	$(CC) $(CFLAGS) ejer2.cpp

ejer3.i :
	$(CC) $(CFLAGS) ejer3.cpp

clean:
	\rm *.o ejer?
