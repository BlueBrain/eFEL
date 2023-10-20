#ifndef EFEL_EXCEPTIONS_H
#define EFEL_EXCEPTIONS_H

#include <stdexcept>

// Define the custom exception class
class EfelAssertionError : public std::runtime_error {
public:
  explicit EfelAssertionError(const std::string& message) : std::runtime_error(message) {}
};

#endif // EFEL_EXCEPTIONS_H
