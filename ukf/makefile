# ===============================
#  Makefile for Generic UKF (C)
# ===============================

# Compiler and flags
CC      := gcc
CFLAGS  := -std=c11 -O2
LDFLAGS := -lm -lc

# Targets
LIB_OBJS := ukf.o
TEST_OBJS := ukf_test.o

# Executable name
TEST_BIN := ukf_test

# Default build
all: $(TEST_BIN)

# -------------------------------
# Build Rules
# -------------------------------
$(TEST_BIN): $(LIB_OBJS) $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

# -------------------------------
# Convenience Targets
# -------------------------------
# Run the demo & print output
run: $(TEST_BIN)
	./$(TEST_BIN)

# Generate plots using Python
plot: $(TEST_BIN)
	./$(TEST_BIN) > ukf_output.txt
	python3 plot_ukf_test.py ukf_output.txt

# Remove object files and binary
clean:
	rm -f *.o $(TEST_BIN) ukf_output.txt

# Full reset (includes backup files)
distclean: clean
	rm -f *~ *.bak

.PHONY: all run plot clean distclean

