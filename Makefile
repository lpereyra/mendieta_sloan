### snapshot options #######
EXTRAS += -DPRINT_XYZ
#EXTRAS += -DLIM_VOL
<<<<<<< HEAD
#EXTRAS += -DLEN_FOF_MERCHAN
=======
EXTRAS += -DLEN_FOF_MERCHAN
>>>>>>> e7d6618c604010219135380dce24974e7baa4b1b

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=8
DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -g
GSLL   := -lgsl -lgslcblas
LIBS   := -lm $(GSLL)

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := grid.o variables.o leesloan.o iden.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf mendieta.o
	rm -rf $(EXEC)
