SUBDIRS = vns smmb_aco
.PHONY: clean
.PHONY: directories

all:
	@for i in $(SUBDIRS); do \
	echo "\033[33m"+"make all in $$i..."+"\033[0m"; \
	(cd $$i; $(MAKE) all); done

clean:
	@for i in $(SUBDIRS); do \
	echo "\033[33m"+"cleaning in $$i..."+"\033[0m"; \
	(cd $$i; $(MAKE) clean); done

install:
	@for i in $(SUBDIRS); do \
	echo "\033[33m"+"Creating directories in $$i..."+"\033[0m"; \
	(cd $$i; $(MAKE) install); done
