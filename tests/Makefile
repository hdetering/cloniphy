CC := g++
CXXFLAGS:=-std=c++11 -W -Wall -pedantic
CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_DEBUG=1

LIB:=-pthread -lboost_program_options -lyaml-cpp -L ../lib #-lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
INC:=-I../include -I../src # place external header files in 'include' directory

TARGET:=../bin/run_tests
SRCDIR:=../src
SRCEXT:=cpp
#SOURCES:=$(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
#SOURCES:=$(filter-out $(SRCDIR)/cloniphy.$(SRCEXT),$(SOURCES))
#OBJECTS:=$(patsubst $(SRCDIR)/%,../debug/%,$(SOURCES:.$(SRCEXT)=.o))
#TEST_SOURCES:=$(shell find . -type f -name '*.$(SRCEXT)')
#TEST_OBJECTS:=$(TEST_SOURCES:.$(SRCEXT)=.o)
SOURCES:=$(shell find . -type f -name '*.$(SRCEXT)')
OBJECTS:=$(shell find ../debug -type f -name "*.o")
OBJECTS:=$(filter-out ../debug/cloniphy.o,$(OBJECTS))
OBJECTS+=$(SOURCES:.$(SRCEXT)=.o)

$(TARGET): $(OBJECTS) #$(TEST_OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB) # doesn't apply later changes to make vars
	#@echo " $(CC) $(OBJECTS) -o $(TARGET) $(LIB)"; $(CC) $(OBJECTS) -o $(TARGET) $(LIB)

#../debug/%.o: $(SRCDIR)/%.$(SRCEXT)
#	@echo " $(CC) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CXXFLAGS) $(INC) -c -o $@ $<

%.o: %.$(SRCEXT)
	@echo " $(CC) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CXXFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) *.o $(TARGET)"; $(RM) *.o $(TARGET)
