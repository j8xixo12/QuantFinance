CXX = g++
INCLUDE = include
SRC = src
OUT = build
EXEC = main

CXXFLAGS = -Iinclude -Wall -g -std=c++17

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
	$(RM) -rf build
	$(RM) $(EXEC)
	cd src/example && $(MAKE) distclean
