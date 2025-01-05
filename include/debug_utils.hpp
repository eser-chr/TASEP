#pragma once

#include <iostream>
#include <stdexcept>

// Debugging and assertion macros
#ifdef DEBUG
    #define ASSERT(condition, message)                                       \
        do {                                                                 \
            if (!(condition)) {                                              \
                std::cerr << "Assertion failed: " << (message)               \
                          << "\nIn file: " << __FILE__                       \
                          << "\nOn line: " << __LINE__ << std::endl;         \
                throw std::runtime_error("Error: " + std::string(message));  \
            }                                                                \
        } while (false)
#else
    #define ASSERT(condition, message)
#endif

// Example debug print macro
#ifdef DEBUG
    #define DEBUG_PRINT(msg) std::cout << "DEBUG: " << (msg) << std::endl
#else
    #define DEBUG_PRINT(msg)
#endif
