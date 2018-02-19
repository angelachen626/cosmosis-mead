include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: libhmcode18.a mead_interface.so

mead_interface.so: libhmcode18.a mead_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ -L. -lhmcode18 $(LDFLAGS) -lcosmosis_fortran -lcosmosis



test: test.f90 libhmcode18.a
	$(FC) $(FFLAGS) -o $@ $< -L. -lhmcode18 $(LDFLAGS)


clean:
	rm -f hmcode18
	rm -f libhmcode18.a
	rm -f hmcode18.o
	rm -f mead_interface.so
	rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod

libhmcode18.a: HMcode18.f90
	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS)
	$(AR) rc $@ HMcode18.o
