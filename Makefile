VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := project

CC := g++-14

ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

ARCH := native

GCC_VER := 14
GCC_PREFIX := /opt/homebrew/Cellar/gcc/$(GCC_VER).2.0_1

CFLAGS := -std=c++17 -O3 -march=$(ARCH) -fopenmp -m64 \
-I$(GCC_PREFIX)/include \
-I$(GCC_PREFIX)/lib/gcc/$(GCC_VER)/gcc/$(shell uname -m)-apple-darwin$(shell uname -r | cut -d. -f1)/$(GCC_VER)/include

COMPILE_COMMAND := $(CC) $(CFLAGS)

BioFVM_OBJECTS := BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o \
BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o 

PhysiCell_core_OBJECTS := PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o \
PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_constants.o PhysiCell_basic_signaling.o \
PhysiCell_signal_behavior.o PhysiCell_rules.o

PhysiCell_module_OBJECTS := PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o \
PhysiCell_pugixml.o PhysiCell_settings.o PhysiCell_geometry.o

PhysiCell_custom_module_OBJECTS := custom.o

pugixml_OBJECTS := pugixml.o

PhysiCell_OBJECTS := $(BioFVM_OBJECTS)  $(pugixml_OBJECTS) $(PhysiCell_core_OBJECTS) $(PhysiCell_module_OBJECTS)
ALL_OBJECTS := $(PhysiCell_OBJECTS) $(PhysiCell_custom_module_OBJECTS)

all: main.cpp $(ALL_OBJECTS)
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp 
	make name

name:
	@echo ""
	@echo "Executable name is" $(PROGRAM_NAME)
	@echo ""

PhysiCell_phenotype.o: ./core/PhysiCell_phenotype.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_phenotype.cpp

PhysiCell_digital_cell_line.o: ./core/PhysiCell_digital_cell_line.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_digital_cell_line.cpp

PhysiCell_cell.o: ./core/PhysiCell_cell.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell.cpp 

PhysiCell_cell_container.o: ./core/PhysiCell_cell_container.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_cell_container.cpp 

PhysiCell_standard_models.o: ./core/PhysiCell_standard_models.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_standard_models.cpp 

PhysiCell_utilities.o: ./core/PhysiCell_utilities.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_utilities.cpp 

PhysiCell_custom.o: ./core/PhysiCell_custom.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_custom.cpp 

PhysiCell_constants.o: ./core/PhysiCell_constants.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_constants.cpp 

PhysiCell_signal_behavior.o: ./core/PhysiCell_signal_behavior.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_signal_behavior.cpp 

PhysiCell_rules.o: ./core/PhysiCell_rules.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_rules.cpp 

BioFVM_vector.o: ./BioFVM/BioFVM_vector.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_vector.cpp 

BioFVM_agent_container.o: ./BioFVM/BioFVM_agent_container.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_agent_container.cpp 

BioFVM_mesh.o: ./BioFVM/BioFVM_mesh.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_mesh.cpp 

BioFVM_microenvironment.o: ./BioFVM/BioFVM_microenvironment.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_microenvironment.cpp 

BioFVM_solvers.o: ./BioFVM/BioFVM_solvers.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_solvers.cpp 

BioFVM_utilities.o: ./BioFVM/BioFVM_utilities.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_utilities.cpp 

BioFVM_basic_agent.o: ./BioFVM/BioFVM_basic_agent.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_basic_agent.cpp 

BioFVM_matlab.o: ./BioFVM/BioFVM_matlab.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_matlab.cpp

BioFVM_MultiCellDS.o: ./BioFVM/BioFVM_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/BioFVM_MultiCellDS.cpp

pugixml.o: ./BioFVM/pugixml.cpp
	$(COMPILE_COMMAND) -c ./BioFVM/pugixml.cpp

PhysiCell_SVG.o: ./modules/PhysiCell_SVG.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_SVG.cpp

PhysiCell_pathology.o: ./modules/PhysiCell_pathology.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pathology.cpp

PhysiCell_MultiCellDS.o: ./modules/PhysiCell_MultiCellDS.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_MultiCellDS.cpp

PhysiCell_various_outputs.o: ./modules/PhysiCell_various_outputs.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_various_outputs.cpp

PhysiCell_pugixml.o: ./modules/PhysiCell_pugixml.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_pugixml.cpp

PhysiCell_settings.o: ./modules/PhysiCell_settings.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_settings.cpp

PhysiCell_basic_signaling.o: ./core/PhysiCell_basic_signaling.cpp
	$(COMPILE_COMMAND) -c ./core/PhysiCell_basic_signaling.cpp 

PhysiCell_geometry.o: ./modules/PhysiCell_geometry.cpp
	$(COMPILE_COMMAND) -c ./modules/PhysiCell_geometry.cpp 

custom.o: ./custom_modules/custom.cpp 
	$(COMPILE_COMMAND) -c ./custom_modules/custom.cpp

clean:
	@rm -f *.o || true
	@rm -f BioFVM/*.o || true
	@rm -f $(PROGRAM_NAME)* || true

.PHONY: all clean name

