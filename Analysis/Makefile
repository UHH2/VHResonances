CC = g++
CFLAGS = -std=c++0x -Wall -g

USERCXXFLAGS := -I${CMSSW_BASE}/src/
USERLDFLAGS  := -L${CMSSW_BASE}/lib/${SCRAM_ARCH} -lHiggsAnalysisCombinedLimit
ROOTLDFLAGS  := -lMinuit -lRooFit -lRooFitCore -L ${ROOTSYS}/lib -lASImage
LFLAGS = -I. -I./include -lm $(USERCXXFLAGS) $(USERLDFLAGS) $(ROOTLDFLAGS)
LIBSMAIN = `root-config --cflags --evelibs --glibs`

INCDIR   = include
SRCDIR   = src
OBJDIR   = obj

SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
SOURCES_UTILS  := $(wildcard $(SRCDIR)/tdr*.cxx)
SOURCES  := $(subst ${SOURCES_UTILS},,${SOURCES})
# INCLUDES := $(SOURCES:$(SRCDIR)/%.cxx=$(INCDIR)/%.hpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
BINS     := $(SOURCES:$(SRCDIR)/%.cxx=%)

OBJECTS_UTILS  := $(SOURCES_UTILS:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
BINS_UTILS     := $(SOURCES_UTILS:$(SRCDIR)/%.cxx=%)


# all: print $(BINS)
all: $(BINS)

print:
	@echo OBJECTS $(OBJECTS)
	@echo OBJECTS_UTILS $(OBJECTS_UTILS)

$(BINS): $(OBJECTS) $(OBJECTS_UTILS)
	@$(CC) -o $@ $(OBJDIR)/$@.o $(OBJECTS_UTILS) $(LFLAGS) $(LIBSMAIN)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	@$(CC) -o $@ -c $^ $(CFLAGS) $(LIBSMAIN) $(LFLAGS)

$(OBJECTS_UTILS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	@$(CC) -o $@ -c $^ $(CFLAGS) $(LIBSMAIN) $(LFLAGS)

clean:
	@rm -f $(OBJECTS)

.PHONY: all clean
