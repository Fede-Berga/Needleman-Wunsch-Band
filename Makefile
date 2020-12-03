main.exe : sw.o
	g++ sw.o -o main.exe

sw.o : sw.c
	g++ -c sw.c -o sw.o

.PHONY: clean

clean:
	rm -rf *.o *.exe