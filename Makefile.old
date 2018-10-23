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

BOOST_FOLDER = /home/storm/Documents/Master2/bitbucket/boost_1_61_0
INC=-I$(BOOST_FOLDER) \
 -I./include \
 -I./include/parsing \
 -I./include/smmb_aco


 CXX=g++ -o
 CXXFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -g $(INC)
 LFLAGS=-std=c++11 -fopenmp -Wall -Wextra -DNDEBUG -I$(BOOST_FOLDER) -lm

 SRCDIR=src
 OBJDIR=obj
 BINDIR=.
 rm=rm -f

 OBJ=$(OBJDIR)/main.o \
 $(OBJDIR)/smmb_aco.o \
 $(OBJDIR)/parsing.o

TARGET=$(BINDIR)/smmb_aco

# all: $(TARGET)

# $(TARGET): $(OBJ)
# 	@g++ $^ -o $@ $(LFLAGS)
# 	@echo "Linking complete."

#---------------------
# main
#---------------------
$(OBJDIR)/main.o: ./$(SRCDIR)/main.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled main.o"

#---------------------
# parsing
#---------------------
$(OBJDIR)/parsing.o: ./$(SRCDIR)/parsing/parsing.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled Parameters_file_parsing.o"

#---------------------
# smmb_aco
#---------------------
$(OBJDIR)/smmb_aco.o: ./$(SRCDIR)/smmb_aco/smmb_aco.cpp
	@$(CXX) $@ -c $< $(CXXFLAGS)
	@echo "Compiled smmb_aco.o"
