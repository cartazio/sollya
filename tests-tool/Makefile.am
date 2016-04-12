TESTFILES = timecheck.sollya check001.sollya check002.sollya check003.sollya check004.sollya check005.sollya check006.sollya check007.sollya check008.sollya check009.sollya check010.sollya check011.sollya check012.sollya check013.sollya check014.sollya check015.sollya check016.sollya check017.sollya check018.sollya check019.sollya check020.sollya check021.sollya check022.sollya check023.sollya check024.sollya check025.sollya check026.sollya check027.sollya check028.sollya check029.sollya check030.sollya

EXTRA_DIST = check.proto $(TESTFILES) $(TESTFILES:.sollya=.sollya.reference)

TESTS = $(TESTFILES:.sollya=.sh)

check_PROGRAMS = force_compilation external_structural_cmp
check_LIBRARIES = libraryexample.a # libraryexample.a will be a shared library actually,
                                   # but we do not know how to explain it to automake.

force_compilation_SOURCES =

MOSTLYCLEANFILES = $(TESTS) libraryexample.a

libraryexample.$(OBJEXT): libraryexample.c
	$(CC) $(CFLAGS) $(AM_CPPFLAGS) -fPIC -Wall -c '$(srcdir)/libraryexample.c'

libraryexample.a: libraryexample.$(OBJEXT)
	$(CC) $(LDFLAGS) $(AM_LDFLAGS) -shared -o libraryexample.a libraryexample.$(OBJEXT) -lmpfi -lmpfr -lgmp

.PHONY: force_compilation$(EXEEXT)

force_compilation$(EXEEXT):
	@echo "Recompiling sollya";cd $(top_builddir); $(MAKE) sollya$(EXEEXT)

external_structural_cmp.o: external_structural_cmp.c
	$(CC) $(CFLAGS) $(AM_CPPFLAGS) -fPIC -Wall -c '$(srcdir)/external_structural_cmp.c'

external_structural_cmp: external_structural_cmp.o
	$(CC) -fPIC -shared -o external_structural_cmp external_structural_cmp.o

$(TESTS): check.proto
	@echo "Building "$@
	@echo "#!/bin/sh" > $@
	@echo "SOLLYA="$(top_builddir)"/sollya" >> $@
	@echo "TESTFILES=\""$(@:.sh=.sollya)"\"" >> $@
	@cat "$(srcdir)/check.proto" >> $@
	@chmod +x $@

