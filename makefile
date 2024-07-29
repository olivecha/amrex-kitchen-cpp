MPICC := mpic++ # This is the main compiler
CC := mpic++
SRCDIR := src
BUILDDIR := build
TARGET := mandoline
  
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CFLAGS :=  -Wall -pedantic -g 

$(TARGET): $(OBJECTS)
		@echo " Linking..." 
			@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
		@mkdir -p $(BUILDDIR)
			@echo " $(MPICC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
		@echo " Cleaning..."; 
			@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
