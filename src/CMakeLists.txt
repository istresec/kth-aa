set(BINARY ${CMAKE_PROJECT_NAME})
add_executable(${BINARY}_kattis main_kattis.cpp TSP.h utility.h blossom5_all_in_one_file.h chokolino.h christofides.h
        local2opt.h clarke_wright.h greedy.h bruteforce.h local2opt_no_knn.h local3opt_no_knn_sequential.h)
add_executable(${BINARY}_test main.cpp TSP.h utility.h blossom5_all_in_one_file.h chokolino.h christofides.h
        local2opt.h clarke_wright.h greedy.h bruteforce.h local2opt_no_knn.h local3opt_no_knn_sequential.h local3opt.h)
