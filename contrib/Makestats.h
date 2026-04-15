#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <iomanip> // std::get_time miatt
#include "test_tools.h"

class Makestats
{
private:
    // Tartományok: [10-19], [20-29], [30-39], [40-49], [50+]
    std::vector<int> thresholds = {10, 20, 30, 40, 50};

    // A csoportosított adatok: Alsó_határ -> Entry-k listája
    std::map<int, std::vector<test_tools::entry>> grouped_data;

    std::time_t stringToTime(const std::string& date_str);

public:
    // Alapértelmezett érték: 2026. január 1. 00:00:00
    Makestats(const std::string& start_date_str = "2026-01-01 00:00:00");
    ~Makestats();

    void make_heu_table();
    void make_heu_chart();
    void make_BnC_table();
    void make_HK_table();
    void make_HK_BnC_table();
};