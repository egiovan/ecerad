FC = gfortran
INCLUDE_DIR = -I../local/include
LIB_DIR = -L../local/lib -lecerad


ecerad_py_interface.so : ecerad_py_interface.f90
	$(FC) -shared -o $@ $< $(INCLUDE_DIR) $(LIB_DIR)


clean:
	rm *.so
