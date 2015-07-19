CC = gcc
NAME = liu_xiuwen
LIBS = -lglut32 -lglu32 -lopengl32
CFLAGS = -O3

LIBS =  -L/usr/X11R6/lib/ -O2 -lglut -lGLU -lGL -lXmu -lXt -lSM -lICE -lXext -lX11 -lXi -lXext -lX11 -lm

lab3: lab3.o SSD_util.o
	$(CC) -o lab3 lab3.o SSD_util.o $(LIBS)
lab3_extra: lab3_extra.o SSD_util.o
	$(CC) -o lab3_extra lab3_extra.o SSD_util.o $(LIBS)
.c.o: 
	$(CC)  $(CFLAGS) -c  $(COPT) $<
tar:
	tar cvfz lab3_$(NAME).tar.gz *.c *.h 
	ls -l lab3_$(NAME).tar.gz
run:
	./lab3 lab3_scene1.ssd
	./lab3 lab3_scene2.ssd
	./lab3 lab3_scene3.ssd
	./lab3 lab3_scene4.ssd

clean:
	rm  *.o

