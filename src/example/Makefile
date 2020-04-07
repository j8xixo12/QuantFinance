CXX = g++
INCLUDE = include
SRC = src
OUT = build
EXEC = 

CXXFLAGS = -I../../include -Wall -g -std=c++17

.PHONY: BIN obj clean distclean

# m

OBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

BINS = $(patsubst $(OUT)/%.o, $(OUT)/%, $(wildcard $(OUT)/*.o))

OBJS := $(addprefix $(OUT)/,$(OBJS))
# deps := $(OBJS:%.o=%.o.d)
# -include $(deps)

$(OUT):
	mkdir -p $(OUT)

$(OUT)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $? -o $@

$(OUT)/%: $(OUT)/%.o
	$(CXX) $(LDFLAGS) $? -o $@

BIN: $(BINS)

all: $(OUT) $(OBJS) BIN 

obj: $(OBJS)

clean:
	${RM} $(OBJS) $(EXEC) $(deps)

distclean: clean
	$(RM) -rf build

