# Makefile for compiling symnmf

CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
TARGET = symnmf
SRC = symnmf.c symnmf_impl.c vector.c matrix.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) -lm

clean:
	rm -f $(TARGET)
