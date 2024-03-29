#===============================================================================
# Makefile for Unix & Linux Systems using a GNU C++ compiler
#===============================================================================

#===============================================================================
# compiler flags
#===============================================================================
# -g --Enable debugging
# -Wall --Turn on all warnings
# -fopenmp to use parallelized computation
# -Wall for all warnings
# -Wextra for extra warnings
# -DNDEBUG for using debugging tools such as dnd or valgrind and profiling of our program if those informations are don't wanted please use -s instead
# -I$(BOOST_FOLDER) include all libraries
# -O1 -O2 -O3  offers three levels of optimization, from least optimized to fastest. The program generated will generally be larger and, above all, the compilation time will increase as you move up levels.

# Add external library path
BOOST_FOLDER = ../../boost_1_61_0
# Give path of .hpp files
INC=-I$(BOOST_FOLDER) \
 -I./include \
 -I./include/parsing \
 -I./include/statistics \
 -I./include/vns

# Save compiler call
 CXX=g++ -o
 #Intial call
 # Flags for less verbose, include parallelized compute, optimization at level 2
 CXXFLAGS=-std=c++11 -fopenmp -s -O2 -g $(INC)
 LFLAGS=-std=c++11 -fopenmp -I$(BOOST_FOLDER) -lm
 # Flags with maximum verbose, include parallelized compute, optimization at level 2, debugger informations
 #CXXFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -O2 -g $(INC)
 #LFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -I$(BOOST_FOLDER) -lm
# Save directory names, store some shell commands
 SRCDIR=src
 OBJDIR=obj
 BINDIR=.
 LOGDIR=log
 rm=rm -f
 mkdir=mkdir -p
#Give .o files to create
OBJ=$(OBJDIR)/main.o \
$(OBJDIR)/parameters_parsing.o \
$(OBJDIR)/file_parsing.o \
$(OBJDIR)/statistics.o \
$(OBJDIR)/contingencies.o \
$(OBJDIR)/vns.o

TARGET=$(BINDIR)/vns.exe

all: $(TARGET)

$(TARGET): $(OBJ)
	@g++ $^ -o $@ $(LFLAGS)
	@echo "\033[32m""Compilation complete, executable produced.""\033[0m"

#===============================================================================
# Main
#===============================================================================
$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled main.o""\033[0m"

#===============================================================================
# VNS
#===============================================================================
$(OBJDIR)/vns.o: ./$(SRCDIR)/vns/vns.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled vns.o""\033[0m"

#===============================================================================
# Parameters_parsing
#===============================================================================
$(OBJDIR)/parameters_parsing.o: ./$(SRCDIR)/parsing/parameters_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled parameters_parsing.o""\033[0m"

#===============================================================================
# file_parsing
#===============================================================================
$(OBJDIR)/file_parsing.o: ./$(SRCDIR)/parsing/file_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled file_parsing.o""\033[0m"

#===============================================================================
# Statistics
#===============================================================================
$(OBJDIR)/statistics.o: ./$(SRCDIR)/statistics/statistics.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled statistics.o""\033[0m"

#===============================================================================
# Statistics
#===============================================================================
$(OBJDIR)/contingencies.o: ./$(SRCDIR)/statistics/contingencies.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled contingencies.o""\033[0m"

#===============================================================================
# PHONEY
#===============================================================================
.PHONEY: clean
clean:
	@$(rm) $(OBJ)
	@echo "\n\033[32m""Cleanup complete.""\033[0m\n"

.PHONEY: remove
remove: clean
	@$(rm) $(TARGET)
	@echo "\n\033[32m""Executable removed.""\033[0m\n"

.PHONEY: install
install:
	@$(mkdir) $(OBJDIR)
	@$(mkdir) $(SRCDIR)
	@$(mkdir) $(LOGDIR)
	@echo "\n\033[32m""Directories created""\033[0m\n"
