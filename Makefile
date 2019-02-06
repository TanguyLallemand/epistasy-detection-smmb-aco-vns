SUBDIRS = vns smmb_aco
.PHONY: clean
.PHONY: directories

all:
	@for i in $(SUBDIRS); do \
	echo "make all in $$i..."; \
	(cd $$i; $(MAKE) all); done

clean:
	@for i in $(SUBDIRS); do \
	echo "cleaning in $$i..."; \
	(cd $$i; $(MAKE) clean); done

directories:
	@for i in $(SUBDIRS); do \
	echo "Creating directories in $$i..."; \
	(cd $$i; $(MAKE) directories); done
