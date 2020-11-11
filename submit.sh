# src/main_kattis.cpp
# src/TSP.h
# src/utility.h
# src/utility.cpp
# src/blossom.h

python2 minifier.py "src/blossom5_all_in_one_file.h" > "tmp/blossom5_all_in_one_file.h"
python2 minifier.py "src/chokolino.h" > "tmp/chokolino.h"
python2 minifier.py "src/christofides.h" > "tmp/christofides.h"
python2 minifier.py "src/main_kattis.cpp" > "tmp/main_kattis.cpp"
python2 minifier.py "src/utility.h" > "tmp/utility.h"
python2 minifier.py "src/bruteforce.h" > "tmp/bruteforce.h"
python2 minifier.py "src/local2opt.h" > "tmp/local2opt.h"
python2 minifier.py "src/local3opt_no_knn_sequential.h" > "tmp/local3opt_no_knn_sequential.h"
python2 minifier.py "src/TSP.h" > "tmp/TSP.h"

python submit.py -p tsp -l C++ tmp/chokolino.h tmp/christofides.h tmp/main_kattis.cpp tmp/utility.h tmp/bruteforce.h tmp/local2opt.h tmp/local3opt_no_knn_sequential.h tmp/TSP.h "tmp/blossom5_all_in_one_file.h"
