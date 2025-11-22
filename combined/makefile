# ===========================================
#  Makefile for UKF + ECM Combined Unit Test
# ===========================================

# Compiler and flags
CC      := gcc
CFLAGS  := -std=c11 -O2
LDFLAGS := -lm -lc

# Sources / objects
SRCS    := ukf.c ecm.c ukf_ecm_test.c
OBJS    := $(SRCS:.c=.o)
BIN     := ukf_ecm_test

# Default target
all: $(BIN)

# -------------------------------------------
# Build rules
# -------------------------------------------
$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Object file rules
ukf.o: ukf.c ukf.h
	$(CC) $(CFLAGS) -c $<

ecm.o: ecm.c ecm.h
	$(CC) $(CFLAGS) -c $<

ukf_ecm_test.o: ukf_ecm_test.c ukf.h ecm.h
	$(CC) $(CFLAGS) -c $<

# -------------------------------------------
# Convenience targets
# -------------------------------------------

# Run the combined test and print CSV to stdout
run: $(BIN)
	./$(BIN)

# Run, save CSV, and (optionally) plot with a Python script
# Adjust plot script name if different.
plot: $(BIN)
	./$(BIN) > ukf_ecm_output.csv
	@if [ -f plot_ukf_ecm.py ]; then \
		echo "Plotting with plot_ukf_ecm.py..."; \
		python3 plot_ukf_ecm.py ukf_ecm_output.csv; \
	else \
		echo "ukf_ecm_output.csv generated (no plot_ukf_ecm.py found)."; \
	fi

# Clean build artifacts
clean:
	rm -f $(OBJS) $(BIN) ukf_ecm_output.csv

distclean: clean
	rm -f *~ *.bak

.PHONY: all run plot clean distclean

