CXX = g++
INCLUDE = include
SRC = src
OUT = build
EXEC = main

CXXFLAGS = -Iinclude -Wall -g -std=c++17
LDFLAGS = -larmadillo

R_HOME := $(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := $(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := $(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := $(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := $(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := $(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CXXFLAGS += $(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL)
LDLIBS = $(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)
LDFLAGS += $(LDLIBS)

.PHONY: run obj clean distclean

all: $(EXEC)

OBJS = $(patsubst $(SRC)/%.cpp, %.o, $(wildcard $(SRC)/*.cpp))

OBJS := $(addprefix $(OUT)/,$(OBJS))
deps := $(OBJS:%.o=%.o.d)
-include $(deps)

$(OUT):
	mkdir -p $(OUT)

$(OUT)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@ 

$(EXEC): $(OUT) $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@ 

example: 
	cd src/example && $(MAKE) all
	cd src/example && $(MAKE) all

run: $(EXEC)
	./$(EXEC)

obj: $(OBJS)

clean:
	${RM} $(OBJS) $(EXEC) $(deps)

distclean: clean
	$(RM) *.dat
	$(RM) *.png
	$(RM) -rf build
	$(RM) $(EXEC)
	cd src/example && $(MAKE) distclean
