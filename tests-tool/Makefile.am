TESTFILES= toto.sollya toto2.sollya
#timecheck.sollya check001.sollya check002.sollya check003.sollya check004.sollya check005.sollya check006.sollya check007.sollya check008.sollya check009.sollya check010.sollya check011.sollya check012.sollya check013.sollya check014.sollya check015.sollya check016.sollya check017.sollya check018.sollya check019.sollya check020.sollya

EXTRA_DIST= check.proto $(TESTFILES) $(TESTFILES:.sollya=.sollya.reference)

TESTS = $(TESTFILES:.sollya=.sh)

.PHONY: $(top_builddir)/sollya

$(top_builddir)/sollya:
	@echo "Recompiling sollya";cd $(top_builddir); $(MAKE) all

$(TESTS): check.proto $(top_builddir)/sollya
	@echo "Building "$@
	@echo "#!/bin/sh" > $@
	@echo "SOLLYA="$(top_builddir)"/sollya" >> $@
	@echo "TESTFILES=\""$(@:.sh=.sollya)"\"" >> $@
	@cat check.proto >> $@
	@chmod +x $@