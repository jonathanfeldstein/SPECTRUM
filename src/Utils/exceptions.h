#ifndef SPECTRUM_EXCEPTIONS_H
#define SPECTRUM_EXCEPTIONS_H

#include <string>
#include <exception>
#include <iostream>
////////////////////////////////////////////////////////
// HYPERGRAPH EXCEPTIONS
////////////////////////////////////////////////////////
struct HyperGraphEmptyException : public std::exception{
    const char* what () const throw (){
        return "Hypergraph is empty";
    }
};

struct HyperGraphSizeException : public std::exception{
    const char* what () const throw (){
        return "Hypergraph is too small (< 3 nodes)";
    }
};

struct HyperGraphConnectedException : public std::exception{
    const char* what () const throw (){
        return "Hypergraph is not connected";
    }
};
struct HyperGraphNotImplementedException : public std::exception{
    const char* what () const throw (){
        return "Path Finding on Hypergraphs not implemented";
    }
};


////////////////////////////////////////////////////////
// DATABASE EXCEPTIONS
////////////////////////////////////////////////////////

struct ArgumentsException : public std::exception{
private:
    std::string message_;
public:
    explicit ArgumentsException(const std::string& message) noexcept;
    virtual ~ArgumentsException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};

struct FileNotFoundException : public std::exception{
private:
    std::string message_;
    std::string filename;
public:
    explicit FileNotFoundException(const std::string& filename) noexcept;
    virtual ~FileNotFoundException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};

struct FileNotOpenedException : public std::exception{
private:
    std::string message_;
public:
    explicit FileNotOpenedException(const std::string& filename) noexcept;
    virtual ~FileNotOpenedException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};


struct DatabaseSyntaxException : public std::exception{
private:
    std::string message_;
public:
    explicit DatabaseSyntaxException(const std::string& line) noexcept;
    virtual ~DatabaseSyntaxException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};
struct InfoSyntaxException : public std::exception{
private:
    std::string message_;
public:
    explicit InfoSyntaxException(const std::string& line) noexcept;
    virtual ~InfoSyntaxException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};
////////////////////////////////////////////////////////
// MAPS EXCEPTIONS
////////////////////////////////////////////////////////
struct KeyNotFoundException : public std::exception{
private:
    std::string message_;
    std::string key;
public:
    explicit KeyNotFoundException(const std::string& key) noexcept;
    virtual ~KeyNotFoundException() = default;
    virtual const char* what() const noexcept override {
        return message_.c_str();
    }
};

#endif //SPECTRUM_EXCEPTIONS_H
