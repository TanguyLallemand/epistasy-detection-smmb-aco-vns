# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Makefile for Unix & Linux Systems #
# using a GNU C++ compiler #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# compiler flags
# -g --Enable debugging
# -Wall --Turn on all warnings
# -D_USE_FIXED_PROTOTYPES_
# --Force the compiler to use the correct headers
# -ansi --Don't use GNU ext; do use ansi standard.

# Add external library path
BOOST_FOLDER = ../../boost_1_61_0
# Give path of .hpp files
INC=-I$(BOOST_FOLDER) \
 -I./include \
 -I./include/parsing \
 -I./include/statistics \
 -I./include/tools \
 -I./include/smmb_aco

# Save compiler call
 CXX=g++ -o
 CXXFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -g $(INC)
 LFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -I$(BOOST_FOLDER) -lm
# Save directory names
 SRCDIR=src
 OBJDIR=obj
 BINDIR=.
 rm=rm -f
#Give .o files to create
OBJ=$(OBJDIR)/main.o \
$(OBJDIR)/file_parsing.o \
$(OBJDIR)/smmb_aco.o \
$(OBJDIR)/parameters_parsing.o \
$(OBJDIR)/tools.o \
$(OBJDIR)/statistics.o

TARGET=$(BINDIR)/smmb_aco

all: $(TARGET)

$(TARGET): $(OBJ)
	@g++ $^ -o $@ $(LFLAGS)
	@echo "Linking complete."

#---------------------
# Main
#---------------------
$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled main.o"

#---------------------
# Parsing
#---------------------
$(OBJDIR)/file_parsing.o: ./$(SRCDIR)/parsing/file_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled file_parsing.o"

$(OBJDIR)/parameters_parsing.o: ./$(SRCDIR)/parsing/parameters_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled parameters_parsing.o"
#---------------------
# Smmb_aco
#---------------------
$(OBJDIR)/smmb_aco.o: ./$(SRCDIR)/smmb_aco/smmb_aco.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled smmb_aco.o"

#---------------------
# Statistics
#---------------------
$(OBJDIR)/statistics.o: ./$(SRCDIR)/statistics/statistics.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled smmb_aco.o"

#---------------------
# Global
#---------------------
$(OBJDIR)/global.o: ./$(SRCDIR)/global.hpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled global.o"

#---------------------
# Tools
#---------------------
$(OBJDIR)/tools.o: ./$(SRCDIR)/tools/tools.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled tools.o"




#---------------------
# PHONEY
#---------------------
.PHONEY: clean
clean:
	@$(rm) $(OBJ)
	@echo "Cleanup complete."

.PHONEY: remove
remove: clean
	@$(rm) $(TARGET)
	@echo "Executable removed."
