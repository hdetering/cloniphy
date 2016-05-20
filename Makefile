default: all

all: debug release

debug:
	$(MAKE) -C debug

release:
	$(MAKE) -C release

test:
	$(MAKE) -C debug
	$(MAKE) -C tests

clean:
	$(MAKE) -C debug clean
	$(MAKE) -C release clean

.PHONY: default all debug release clean
