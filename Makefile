CXXFLAGS=-fopenmp -O0

numericalIntegration : numericalIntegration.o
	$(CXX) $(CXXFLAGS) numericalIntegration.o -onumericalIntegration

clean:
	rm -f numericalIntegration *~ *.o
