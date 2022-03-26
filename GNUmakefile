.PHONY: test np

np: 
	$(MAKE) -C $@

test: np/np.so
	g++ -o test_npcpp test_npcpp.cpp np/np.so 

all: np test
