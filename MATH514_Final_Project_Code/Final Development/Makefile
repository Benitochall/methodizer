include ../Make.inc

CFLAGS-add += -I../src -DNDEBUG -w
LDFLAGS-add += -L../src -llsoda -lm

first: textUI
second: ode

textUI.o: textUI.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -c $< -o $@
textUI: textUI.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDFLAGS-add)
ode.o: ode.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) $(fPIC) -c $< -o $@
ode: ode.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDFLAGS-add)
clean:
	rm -f *.o ode textUI ysin1y ode
