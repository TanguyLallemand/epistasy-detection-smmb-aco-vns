#BOOST_FOLDER=/ccipl/logiciels/boost/include/
#BOOST_FOLDER=/home/clem/Softwares/boost_1_63_0/
BOOST_FOLDER=/home/c-sinoquet/downloads/boost_1_61_0/
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

INC=-I$(BOOST_FOLDER) \
 -I./include \
 -I./include/parsing \
 -I./include/services \
 -I./include/structures \
 -I./include/statistics \
 -I./include/information_theory \
 -I./include/usecases \
 -I./include/markov

CXX=g++ -o
CXXFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -g $(INC)
LFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -I$(BOOST_FOLDER) -lm

SRCDIR=src
OBJDIR=obj
BINDIR=.
rm=rm -f

OBJ=$(OBJDIR)/main.o \
$(OBJDIR)/Parameters_file_parsing.o \
$(OBJDIR)/Services.o \
$(OBJDIR)/Contingency.o \
$(OBJDIR)/G2_test_indep.o \
$(OBJDIR)/G2_conditional_test_indep.o \
$(OBJDIR)/Mutual_information.o \
$(OBJDIR)/Permutations.o \
$(OBJDIR)/Permutations_adapt.o \
$(OBJDIR)/Smmb_aco.o \
$(OBJDIR)/Smmb_aco_usecase.o

TARGET=$(BINDIR)/SMMB_ACO

all: $(TARGET)

$(TARGET): $(OBJ)
	@g++ $^ -o $@ $(LFLAGS)
	@echo "Linking complete."

$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled main.o"

#---------------------
# Parsing
#---------------------
$(OBJDIR)/Parameters_file_parsing.o: ./$(SRCDIR)/parsing/Parameters_file_parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Parameters_file_parsing.o"

#---------------------
# Statistics
#---------------------
$(OBJDIR)/Contingency.o: ./$(SRCDIR)/statistics/Contingency.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Contingency.o"

$(OBJDIR)/G2_conditional_test_indep.o: ./$(SRCDIR)/statistics/G2_conditional_test_indep.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled G2_conditional_test_indep.o"

$(OBJDIR)/G2_test_indep.o: ./$(SRCDIR)/statistics/G2_test_indep.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled G2_test_indep.o"

$(OBJDIR)/Permutations.o: ./$(SRCDIR)/statistics/Permutations.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Permutations.o"

$(OBJDIR)/Permutations_adapt.o: ./$(SRCDIR)/statistics/Permutations_adapt.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Permutations_adapt.o"

#---------------------
# Services
#---------------------
$(OBJDIR)/Services.o: ./$(SRCDIR)/services/Services.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Services.o"

#---------------------
# Markov
#---------------------
$(OBJDIR)/Smmb_aco.o: ./$(SRCDIR)/markov/Smmb_ACO.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Smmb_aco.o"

#---------------------
# Information theory
#---------------------
$(OBJDIR)/Mutual_information.o: ./$(SRCDIR)/information_theory/Mutual_information.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Mutual_information.o"

#---------------------
# Usecases
#---------------------
$(OBJDIR)/Smmb_aco_usecase.o: ./$(SRCDIR)/usecases/Smmb_aco_usecase.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Smmb_aco_usecase.o"

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
