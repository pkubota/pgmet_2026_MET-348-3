FC = gfortran
FFLAGS = -O2 -Wall -Wextra -std=f2008 -Wno-integer-division

all: dfi_advdif saida_py/campos_t0.dat saida_f90/campos_t0.dat
	python3 plot_resultados.py

dfi_advdif: dfi_advdif.f90
	$(FC) $(FFLAGS) -o $@ $<

saida_py/campos_t0.dat: dfi_advdif.py
	python3 dfi_advdif.py --outdir saida_py

saida_f90/campos_t0.dat: dfi_advdif
	./dfi_advdif saida_f90

clean:
	rm -rf dfi_advdif saida_py saida_f90 resultados_dfi_advdif.png
