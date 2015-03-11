CC = mpicc
CFLAGS = -Wall -o2
TARGET = primegap
all: $(TARGET)

$(TARGET): primegap.c
	$(CC) $(CFLAGS) -o $(TARGET) primegap.c -lgmp
clean:
	$(RM) $(TARGET)
