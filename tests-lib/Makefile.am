EXTRA_DIST=check00.c

#TESTS_ENVIRONMENT=$(SHELL) 

TESTS= check00

.PHONY: $(top_builddir)/sollya

$(top_builddir)/sollya:
	cd $(top_builddir); make all

check00.o: check00.c
	$(CC) check00.c $(CFLAGS) -c

check00: check00.o
	$(CC) check00.o $(LDLAGS) -o check00