CC=gcc
CC+=-DDEBUG -g
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.C EventReader.C EventReaderDict.C 
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=AnalysisLSAlphas

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) -lstdc++

CFLAGS += -I../ 

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.C.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o; rm $(EXECUTABLE)
