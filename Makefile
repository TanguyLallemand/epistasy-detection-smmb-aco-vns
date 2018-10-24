SUBDIRS = vns smmb_aco
all:
	@for i in $(SUBDIRS); do \
	echo "make all in $$i..."; \
	(cd $$i; $(MAKE) all); done
