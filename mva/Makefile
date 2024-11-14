CXX = g++
CPPFLAGS = -g -Wall -O2 -lstdc++fs
CPPFLAGS += -I include/

ROOTFLAGS = $(shell root-config --libs --cflags)
ROOTFLAGS += -lMLP -lMinuit -lTreePlayer -lTMVA -lXMLIO  -lMLP -lm -lTMVAUtils

MAKEFLAGS := --jobs=8

SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,build/%.o,$(SRCS))

TARGET = runMVA

BUILD_DIR = build/

all: oogabooga

oogabooga: $(OBJS)
	$(CXX) $(CPPFLAGS) $(MYINCLUDES) $(ROOTFLAGS) $(OBJS) -o bin/$(TARGET)

build/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(MYINCLUDES) $(ROOTFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.DELETE_ON_ERROR:

clean:
	$(RM) $(OBJS) bin/$(TARGET)

clearplots:
	$(RM)

.PHONY: all oogabooga clean clearplots
