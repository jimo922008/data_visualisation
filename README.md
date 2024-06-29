## Installation
(1) Download files from git repository.

https://gitlab.developers.cam.ac.uk/phy/data-intensive-science-mphil/projects/mj425.git

(2) Generate makefile: by running the command:
```bash
cmake -S . -B build
```

(3) Compile the code: by running the command 
```bash
cmake --build build
```

(4) Run the code: by running the command:
```bash
mpirun -np {number of ranks} ./bin/sheap_mpi {input file} 
```
For example, to run the example LJ13 data set, you can run the command:
```bash
mpirun -np 4 ./bin/sheap_mpi test/LJ13.vec
```
Alternatively, you can also run the command if you would like the code with single rank. 
```bash
./bin/sheap test/LJ13.vec
```

(5) Clean the code: by running the command 
```bash
cmake --build build --target clean
```

Other functions include: 

(6) Generate documentation: by running the command 
```bash
doxygen build/Doxyfile.doc
```

(7) Run tests: by running the command 
```bash
bin/test_optimisaton
```

The coursework has utlised Copilot to generate the some of the functions and tests, all the documentations and the readme file. 