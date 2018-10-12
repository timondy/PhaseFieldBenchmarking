all:	REALSTRIPPER

REALSTRIPPER: realtimeediter.c
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) $(CMD) *.o *.i *.mod *~

