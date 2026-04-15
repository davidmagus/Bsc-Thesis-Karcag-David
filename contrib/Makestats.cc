#include "Makestats.h"
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <algorithm>
#include <climits>

std::time_t Makestats::stringToTime(const std::string &date_str)
{
    struct tm tm = {};
    std::istringstream ss(date_str);
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    return std::mktime(&tm);
}

Makestats::Makestats(const std::string &start_date_str)
{
    std::time_t start_limit = stringToTime(start_date_str);
    std::ifstream file("results.txt");
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Hiba: results.txt nem nyithato meg!" << std::endl;
        return;
    }

    // Fejléc átugrása
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        std::istringstream ss(line);
        std::string d, t;
        ss >> d >> t; // Dátum és idő beolvasása (pl. 2026-04-13 12:31:22)

        if (stringToTime(d + " " + t) < start_limit)
            continue;

        int n_val, seed_val;
        ss >> n_val >> seed_val;

        // Entry létrehozása és adatok betöltése a fájl sorrendjében
        test_tools::entry e(n_val, seed_val);

        ss >> e.NN_val >> e.NN_time >> e.NN_percentage;
        ss >> e.RNN_val >> e.RNN_time >> e.RNN_percentage;
        ss >> e.Grdy_val >> e.Grdy_time >> e.Grdy_percentage;
        ss >> e.Mininsert_val >> e.Mininsert_time >> e.Mininsert_percentage;
        ss >> e.Maxinsert_val >> e.Maxinsert_time >> e.Maxinsert_percentage;
        ss >> e.Randinsert_val >> e.Randinsert_time >> e.Randinsert_percentage;

        std::string bnc_s, hk_s;
        ss >> bnc_s >> e.BnC_val >> e.BnC_Time;
        ss >> hk_s >> e.HK_val >> e.HK_Time;

        e.BnC_solved = (bnc_s == "YES");
        e.HK_solved = (hk_s == "YES");

        // Besorolás a megfelelő 'bucket'-be
        int bucket = 0;
        for (int th : thresholds)
        {
            if (n_val >= th)
                bucket = th;
            else
                break;
        }

        if (bucket > 0)
        {
            grouped_data[bucket].push_back(e);
        }
    }
    file.close();
}

Makestats::~Makestats() {}

// Megjegyzés: A táblázat- és ábrakészítő függvényeket később implementáljuk.

void Makestats::make_heu_table()
{
    std::string html_name = "heu_table.html";
    std::ofstream html(html_name);

    html << "<html><head><style>"
         << "table { border-collapse: collapse; font-family: sans-serif; width: 100%; }"
         << "th, td { border: 1px solid #ccc; padding: 8px; text-align: center; }"
         << "th { background-color: #f2f2f2; }"
         << ".header { background-color: #d1d1d1; font-weight: bold; }"
         << "</style></head><body>"
         << "<h2> Heuristics accuracy (% of OPT)</h2>"
         << "<table>"
         << "<tr><th rowspan='2'>n</th><th colspan='2'>NN</th><th colspan='2'>RNN</th>"
         << "<th colspan='2'>Greedy</th><th colspan='2'>MaxInsert</th>"
         << "<th colspan='2'>MinInsert</th><th colspan='2'>RandInsert</th></tr>"
         << "<tr><th>Avg</th><th>Worst</th><th>Avg</th><th>Worst</th>"
         << "<th>Avg</th><th>Worst</th><th>Avg</th><th>Worst</th>"
         << "<th>Avg</th><th>Worst</th><th>Avg</th><th>Worst</th></tr>";

    // Adatok feldolgozása bucket-enként
    for (auto const &[threshold, entries] : grouped_data)
    {
        if (entries.empty())
            continue;

        double sumNN = 0, worstNN = 0;
        double sumRNN = 0, worstRNN = 0;
        double sumGr = 0, worstGr = 0;
        double sumMaxI = 0, worstMaxI = 0;
        double sumMinI = 0, worstMinI = 0;
        double sumRandI = 0, worstRandI = 0;

        for (const auto &e : entries)
        {
            sumNN += e.NN_percentage;
            worstNN = std::max(worstNN, e.NN_percentage);
            sumRNN += e.RNN_percentage;
            worstRNN = std::max(worstRNN, e.RNN_percentage);
            sumGr += e.Grdy_percentage;
            worstGr = std::max(worstGr, e.Grdy_percentage);
            sumMaxI += e.Maxinsert_percentage;
            worstMaxI = std::max(worstMaxI, e.Maxinsert_percentage);
            sumMinI += e.Mininsert_percentage;
            worstMinI = std::max(worstMinI, e.Mininsert_percentage);
            sumRandI += e.Randinsert_percentage;
            worstRandI = std::max(worstRandI, e.Randinsert_percentage);
        }

        size_t count = entries.size();

        // Sor kiírása
        html << "<tr>"
             << "<td class='header'>" << threshold << "-" << threshold + 10 << "</td>"
             << "<td>" << std::fixed << std::setprecision(2) << (sumNN / count) * 100 << "%</td><td>" << worstNN * 100 << "%</td>"
             << "<td>" << (sumRNN / count) * 100 << "%</td><td>" << worstRNN * 100 << "%</td>"
             << "<td>" << (sumGr / count) * 100 << "%</td><td>" << worstGr * 100 << "%</td>"
             << "<td>" << (sumMaxI / count) * 100 << "%</td><td>" << worstMaxI * 100 << "%</td>"
             << "<td>" << (sumMinI / count) * 100 << "%</td><td>" << worstMinI * 100 << "%</td>"
             << "<td>" << (sumRandI / count) * 100 << "%</td><td>" << worstRandI * 100 << "%</td>"
             << "</tr>";
    }

    html << "</table></body></html>";
    html.close();

    // PNG-vé alakítás (Ubuntu alatt az ImageMagick 'import' vagy 'wkhtmltoimage' ajánlott)
    // Ha nincs fent: sudo apt install wkhtmltopdf
    std::cout << "HTML tabla elkeszult: heu_table.html" << std::endl;

    // Rendszerhívás a képkonverzióhoz (opcionális, ha telepítve van az eszköz)
    // system("wkhtmltoimage heu_table.html heu_table.png");
}
void Makestats::make_heu_chart()
{
    // 1. Adatok előkészítése gnuplot számára
    std::ofstream avg_file("heu_avg.dat");
    std::ofstream worst_file("heu_worst.dat");

    avg_file << "# n NN RNN Greedy MaxI MinI RandI" << std::endl;
    worst_file << "# n NN RNN Greedy MaxI MinI RandI" << std::endl;

    for (auto const &[threshold, entries] : grouped_data)
    {
        if (entries.empty())
            continue;

        double sNN = 0, sRNN = 0, sGr = 0, sMax = 0, sMin = 0, sRand = 0;
        double wNN = 0, wRNN = 0, wGr = 0, wMax = 0, wMin = 0, wRand = 0;

        for (const auto &e : entries)
        {
            sNN += e.NN_percentage;
            wNN = std::max(wNN, e.NN_percentage);
            sRNN += e.RNN_percentage;
            wRNN = std::max(wRNN, e.RNN_percentage);
            sGr += e.Grdy_percentage;
            wGr = std::max(wGr, e.Grdy_percentage);
            sMax += e.Maxinsert_percentage;
            wMax = std::max(wMax, e.Maxinsert_percentage);
            sMin += e.Mininsert_percentage;
            wMin = std::max(wMin, e.Mininsert_percentage);
            sRand += e.Randinsert_percentage;
            wRand = std::max(wRand, e.Randinsert_percentage);
        }

        size_t c = entries.size();
        int mid_n = threshold + 5; // A tartomány közepe a tengelyen

        avg_file << mid_n << " " << (sNN / c) * 100 << " " << (sRNN / c) * 100 << " " << (sGr / c) * 100
                 << " " << (sMax / c) * 100 << " " << (sMin / c) * 100 << " " << (sRand / c) * 100 << std::endl;

        worst_file << mid_n << " " << wNN * 100 << " " << wRNN * 100 << " " << wGr * 100
                   << " " << wMax * 100 << " " << wMin * 100 << " " << wRand * 100 << std::endl;
    }
    avg_file.close();
    worst_file.close();

    std::ofstream gp("heu_plots.gp");
    gp << "set terminal pngcairo size 800,500 font 'sans,10'" << std::endl;
    gp << "set xlabel 'n (Problem Range)'" << std::endl;
    gp << "set ylabel 'Performance (% of Optimum)'" << std::endl;
    gp << "set grid" << std::endl;
    gp << "set key outside right center" << std::endl;
    gp << "set xtics 10" << std::endl;

    // Átlag ábra - Itt heu_avg.dat kell
    gp << "set output 'heu_avg.png'" << std::endl;
    gp << "set title 'Average Performance (Grouped by 10s)'" << std::endl;
    gp << "plot 'heu_avg.dat' u 1:2 w lp pt 7 t 'NN', '' u 1:3 w lp pt 7 t 'RNN', "
       << "'' u 1:4 w lp pt 7 t 'Greedy', '' u 1:5 w lp pt 7 t 'MaxI', "
       << "'' u 1:6 w lp pt 7 t 'MinI', '' u 1:7 w lp pt 7 t 'RandI'" << std::endl;

    // Worst ábra - JAVÍTVA: heu_worst.dat kell, nem worst_avg.dat!
    gp << "set output 'heu_worst.png'" << std::endl;
    gp << "set title 'Worst Case Performance (Grouped by 10s)'" << std::endl;
    gp << "plot 'heu_worst.dat' u 1:2 w lp pt 7 t 'NN', '' u 1:3 w lp pt 7 t 'RNN', "
       << "'' u 1:4 w lp pt 7 t 'Greedy', '' u 1:5 w lp pt 7 t 'MaxI', "
       << "'' u 1:6 w lp pt 7 t 'MinI', '' u 1:7 w lp pt 7 t 'RandI'" << std::endl;
    gp.close();

    system("gnuplot heu_plots.gp");

    // Gnuplot futtatása
    system("gnuplot heu_plots.gp");

    // 3. HTML kiegészítése a képekkel
    std::ofstream html("heu_table.html", std::ios::app);
    html << "<br><hr><h2>Visualization</h2>"
         << "<div style='display: flex; flex-direction: column; align-items: center;'>"
         << "<img src='heu_avg.png' style='margin-bottom: 20px;'>"
         << "<img src='heu_worst.png'>"
         << "</div></body></html>";
    html.close();

    std::cout << "Abrak elmentve es hozzaadva a heu_table.html-hez." << std::endl;
}
void Makestats::make_BnC_table() {}
void Makestats::make_HK_table() {}
void Makestats::make_HK_BnC_table() {}