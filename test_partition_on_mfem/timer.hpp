#ifndef TIMER
#define TIMER

#include <iostream>
#include <chrono>
#include <iomanip>
class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::string timer_name;
    
public:
    Timer(const std::string& name = "Timer") : timer_name(name) {
        start();
    }
    
    ~Timer() {
    }
    
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    void stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        
        std::cout << std::fixed << std::setprecision(3);
        if (duration.count() < 1000) {
            std::cout << timer_name << " took: " << duration.count() << " µs" << std::endl;
        } else if (duration.count() < 1000000) {
            std::cout << timer_name << " took: " << duration.count() / 1000.0 << " ms" << std::endl;
        } else {
            std::cout << timer_name << " took: " << duration.count() / 1000000.0 << " s" << std::endl;
        }
    }
    
    template<typename Duration = std::chrono::milliseconds>
    long long elapsed() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<Duration>(end_time - start_time).count();
    }
    
    void reset() {
        start_time = std::chrono::high_resolution_clock::now();
    }
};

#endif