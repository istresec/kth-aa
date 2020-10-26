#include "gtest/gtest.h"
#include "TSP.h"

#include <vector>

using namespace std;

TEST(naiveTest, test1) {
    //arrange
    //act
    //assert
    vector<int> expected_result = {0, 8, 5, 4, 3, 9, 6, 2, 1, 7};
    vector<pair<double, double>> cities = {
            make_pair(95.0129, 61.5432),
            make_pair(23.1139, 79.1937),
            make_pair(60.6843, 92.1813),
            make_pair(48.5982, 73.8207),
            make_pair(89.1299, 17.6266),
            make_pair(76.2097, 40.5706),
            make_pair(45.6468, 93.5470),
            make_pair(1.8504, 91.6904),
            make_pair(82.1407, 41.0270),
            make_pair(44.4703, 89.3650)
    };

    EXPECT_EQ (TSP::travel_naive(cities), expected_result);
    EXPECT_GE (TSP::travel_naive(cities), expected_result);
}