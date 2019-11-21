CC := g++
SRCDIR := src

BUILDDIR := build
TARGET_TOMO := bin/seismic_tomography

LIBSRC_FMM := $(SRCDIR)/fmmiolib.cpp $(SRCDIR)/fmmlib.cpp $(SRCDIR)/raylib.cpp $(SRCDIR)/velocitylib.cpp $(SRCDIR)/numtools.cpp

SOURCES_TOMO := $(LIBSRC_FMM) $(SRCDIR)/tomolib.cpp $(SRCDIR)/mcmc.cpp
OBJECTS_TOMO := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES_TOMO:.cpp=.o))

CFLAGS := -std=c++11 -O3 
LIB := -larmadillo
INC := -I include

all: $(TARGET_TOMO) 

$(TARGET_TOMO): $(OBJECTS_TOMO)
	@echo " Linking..."
	@echo " $(CC) $^ -o $@ $(LIB)"; $(CC) $^ -o $@ $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET_TOMO) "; $(RM) -r $(BUILDDIR) $(TARGET_TOMO)

.PHONY: clean
