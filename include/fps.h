#pragma once

#include <chrono>
#include <climits>
#include <cstdint>


inline static constexpr auto FPS_UPDATE_TIME = 1000;

class fps_calculator
{
    double reported_fps = 0;
    double reported_min_fps = 0;
    double reported_max_fps = 0;

    int64_t prev_report_tick = 0;
    int64_t prev_tick = 0;
    int64_t frame_count = 0;
    int64_t max_ms = INT_MIN;
    int64_t min_ms = INT_MAX;

public:
    void report_frame()
    {
        int64_t current_tick = duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

        // probably not the unix epoch right now...
        if (prev_tick == 0)
        {
            prev_report_tick = prev_tick = current_tick;
            return;
        }

        // update everything
        frame_count++;
        max_ms = std::max(max_ms, current_tick - prev_tick);
        min_ms = std::min(max_ms, current_tick - prev_tick);
        prev_tick = current_tick;

        // update the reported information
        if (current_tick - prev_report_tick > FPS_UPDATE_TIME)
        {
            reported_fps = frame_count * FPS_UPDATE_TIME / double(current_tick - prev_report_tick);

            // 1/(ms/frame) = frame/ms
            // frame/ms * 1000ms/s = frame/s
            reported_min_fps = FPS_UPDATE_TIME / double(min_ms);
            reported_max_fps = FPS_UPDATE_TIME / double(max_ms);
            prev_report_tick = current_tick;
            max_ms = INT_MIN;
            min_ms = INT_MAX;
            frame_count = 0;
        }
    }

    [[nodiscard]] constexpr auto get_fps() const -> double { return reported_fps; }
    [[nodiscard]] constexpr auto get_min_fps() const -> double { return reported_min_fps; }
    [[nodiscard]] constexpr auto get_max_fps() const -> double { return reported_max_fps; }
};


