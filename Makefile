#all:
#	gcc -std=gnu99 *.c *.h -I. -o gmres  -lm -llapacke -lblas -g -pg -Wall

#clean:
#	rm gmres gmon.out


CC=gcc -std=gnu99
CFLAGS=-g -pg -lm -Wall

OBJS=main.o gmres.o functions.o matrix.o

all: gmres

gmres: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)
		
main.o: main.c
	$(CC) -c $^

gmres.o: gmres.c
	$(CC) -c $< $(CFLAGS)
	
functions.o: functions.c
	$(CC) -c $< $(CFLAGS)

matrix.o: matrix.c 
	$(CC) -c $<
	
clean:
	rm -rf $(OBJS) gmres gmon.out

