xsec: currents_opt_v1.o mathtool.o xsec_fact_new.o main_xsec_new.o \
  nform.o
	gfortran -o $@ $^

%.o: %.f90
	gfortran -o $@ -c $^

.PHONY: clean

clean:
	$(RM) *.mod *.o xsec
