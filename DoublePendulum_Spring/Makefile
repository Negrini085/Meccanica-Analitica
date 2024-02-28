LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

ex: main.cpp EquazioniDifferenziali.h FunzioniVettoriali.h OperazioniVector.h
	g++ main.cpp -o ex ${INCS} ${LIBS}

clean:
	rm ex
