SHELL = /bin/sh
RM = rm -f
LN = ln -fs

DUNE = dune

EXEC = index_calculus
TARGET = calc_index.bc

.PHONY: all $(TARGET) $(EXEC) clean

all: $(TARGET) $(EXEC)

$(TARGET):
	$(DUNE) build $(addprefix src/,$@)

$(EXEC): $(TARGET)
	$(LN) _build/default/src/$< $(EXEC)

clean:
	$(DUNE) clean
	$(RM) $(EXEC)
