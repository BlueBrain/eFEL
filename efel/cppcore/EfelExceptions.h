#ifndef EFEL_EXCEPTIONS_H
#define EFEL_EXCEPTIONS_H

#include <stdexcept>

// Define the custom exception class
// Raise this error to mark the function failed when a required condition is not met 
class EfelAssertionError : public std::runtime_error {
public:
  explicit EfelAssertionError(const std::string& message) : std::runtime_error(message) {}
};

class FeatureComputationError : public std::runtime_error {
public:
  explicit FeatureComputationError(const std::string& message) 
    : std::runtime_error("An error occurred while computing the feature, feature is not found. " + message) {}
};

class EmptyFeatureError : public std::runtime_error {
public:
  explicit EmptyFeatureError(const std::string& message) 
    : std::runtime_error("The feature is found in the feature dictionary but it is empty. " + message) {}
};

#endif // EFEL_EXCEPTIONS_H
