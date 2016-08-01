#makefile

PROGRAM = main.exe
OBJS = main.o ga.o rnd.o kpga.o
CC = gcc

$(PROGRAM): $(OBJS)
	$(CC) -o $(PROGRAM) $(OBJS)

main.o: main.c
	$(CC) -c $<

ga.o: ga.c
	$(CC) -c $<

rnd.o: rnd.c
	$(CC) -c $<

kpga.o: kpga.c
	$(CC) -c $<

clean:
	rm -f $(PROGRAM) $(OBJS)

