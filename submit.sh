# src/main_kattis.cpp
# src/TSP.h
# src/utility.h
# src/utility.cpp
# src/blossom.h

python2 minifier.py "src/blossom.h" > "tmp/blossom.h"

python2 minifier.py "src/main_kattis.cpp" > "tmp/main_kattis.cpp"
# cp "src/main_kattis.cpp" "tmp/main_kattis.cpp"

python2 minifier.py "src/TSP.h" > "tmp/TSP.h"
# cp "src/TSP.h" "tmp/TSP.h"

python2 minifier.py "src/utility.h" > "tmp/utility.h"
# cp "src/utility.h" "tmp/utility.h"

python2 minifier.py "src/utility.cpp" > "tmp/utility.cpp"
# cp "src/utility.cpp" "tmp/utility.cpp"


python submit.py -p tsp -l C++ tmp/main_kattis.cpp tmp/TSP.h tmp/utility.h tmp/utility.cpp tmp/blossom.h
