# =====================================
#  Makefile for ECM Battery Model Test
# =====================================

# Compiler and flags
CC      := gcc
CFLAGS  := -std=c11 -Wall -Wextra -O2
LDFLAGS := -lm -lc

# Sources / objects
ECM_SRCS  := ecm.c ecm_test.c
ECM_OBJS  := $(ECM_SRCS:.c=.o)
ECM_BIN   := ecm_test

# Default target
all: $(ECM_BIN)

# -------------------------------------
# Build rules
# -------------------------------------
$(ECM_BIN): $(ECM_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

# ecm_test.c does not have its own header dependency
ecm_test.o: ecm_test.c ecm.h
	$(CC) $(CFLAGS) -c $<

# -------------------------------------
# Convenience targets
# -------------------------------------

# Run the ECM unit test and print CSV output
run: $(ECM_BIN)
	./$(ECM_BIN)

# Run test, save CSV, and plot via Python
plot: $(ECM_BIN)
	./$(ECM_BIN) > ecm_output.csv
	python3 plot_ecm.py ecm_output.csv

# Remove build artifacts
clean:
	rm -f $(ECM_OBJS) $(ECM_BIN) ecm_output.csv

distclean: clean
	rm -f *~ *.bak

.PHONY: all run plot clean distclean

