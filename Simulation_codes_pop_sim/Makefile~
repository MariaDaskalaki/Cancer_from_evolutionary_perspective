SOURCES =  test_maria_test.c
OBJECTS = $(SOURCES:%.c=%.o)

TARGET = test_maria_test

all: $(TARGET)

%.o: %.c
	#gcc -openmp -pedantic -Wall -Werror  -g $< -c -o $@
	gcc -openmp -pedantic -w  -g $< -c -o $@ ### to -Wall -Werror briskei finds all the warnings and converts all the warnings into errors in order to not compile.

$(TARGET): $(OBJECTS)
	gcc $(OBJECTS) -o $@ -lm -lgsl -lgslcblas

clean:
	rm -f $(TARGET)
	rm -f $(OBJECTS)
	rm -f *.png
