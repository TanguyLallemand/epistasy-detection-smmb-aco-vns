#===============================================================================
# Makefile for Unix & Linux Systems using a GNU C++ compiler #
#===============================================================================

#===============================================================================
# Compiler flags
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
 -I./include/tools \
 -I./include/smmb_aco

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
$(OBJDIR)/file_parsing.o \
$(OBJDIR)/smmb_aco.o \
$(OBJDIR)/parameters_parsing.o \
$(OBJDIR)/tools.o \
$(OBJDIR)/contingencies.o \
$(OBJDIR)/statistics.o



TARGET=$(BINDIR)/smmb_aco.exe

all:$(TARGET)

$(TARGET):$(OBJ)
	@g++ $^ -o $@ $(LFLAGS)
	@echo "\033[32m""Compilation complete, executable produced.""\033[0m"

#===============================================================================
# Main
#===============================================================================
$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled main.o""\033[0m"

#===============================================================================
# Parsing
#===============================================================================
$(OBJDIR)/file_parsing.o: ./$(SRCDIR)/parsing/file_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled file_parsing.o""\033[0m"

$(OBJDIR)/parameters_parsing.o: ./$(SRCDIR)/parsing/parameters_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled parameters_parsing.o""\033[0m"
#===============================================================================
# Smmb_aco
#===============================================================================
$(OBJDIR)/smmb_aco.o: ./$(SRCDIR)/smmb_aco/smmb_aco.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled smmb_aco.o""\033[0m"

#===============================================================================
# Statistics
#===============================================================================
$(OBJDIR)/statistics.o: ./$(SRCDIR)/statistics/statistics.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled statistics.o""\033[0m"

$(OBJDIR)/contingencies.o: ./$(SRCDIR)/statistics/contingencies.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled contingencies.o""\033[0m"

#===============================================================================
# Global
#===============================================================================
$(OBJDIR)/global.o: ./$(SRCDIR)/global.hpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled global.o""\033[0m"

#===============================================================================
# Tools
#===============================================================================
$(OBJDIR)/tools.o: ./$(SRCDIR)/tools/tools.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled tools.o""\033[0m"

#===============================================================================
# Output_writing
#===============================================================================
$(OBJDIR)/output_writing.o: ./$(SRCDIR)/output_writing/output_writing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "\033[32m""Compiled output_writing.o""\033[0m"



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
