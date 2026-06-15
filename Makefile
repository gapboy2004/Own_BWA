CC      = gcc
CFLAGS  = -Wall -Wextra -O2 -g -std=c11
SRCS    = bwt_build.c fm_index.c bwt_smem.c
OBJS    = $(SRCS:.c=.o)

.PHONY: all test clean

all: libbwt.a

libbwt.a: $(OBJS)
	ar rcs $@ $^

test: test_bwt
	./test_bwt

test_bwt: test_bwt.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c bwt.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) libbwt.a test_bwt