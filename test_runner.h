#pragma once

#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;

template <class T>
ostream& operator << (ostream& os, const array<T,3>& s){
	  os << "{";
	  bool first = true;
	  for (const auto& x : s) {
	    if (!first) {
	      os << ", ";
	    }
	    first = false;
	    os << x;
	  }
	  return os << "}";
}

template <class T>
ostream& operator << (ostream& os, const vector<T>& s) {
  os << "{";
  bool first = true;
  for (const auto& x : s) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << x;
  }
  return os << "}";
}

template <class T>
ostream& operator << (ostream& os, const set<T>& s) {
  os << "{";
  bool first = true;
  for (const auto& x : s) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << x;
  }
  return os << "}";
}

template <class K, class V>
ostream& operator << (ostream& os, const map<K, V>& m) {
  os << "{";
  bool first = true;
  for (const auto& kv : m) {
    if (!first) {
      os << ", ";
    }
    first = false;
    os << kv.first << ": " << kv.second;
  }
  return os << "}";
}

template<class T, class U>
void AssertEqual(const T& t, const U& u, const string& hint = {}) {
  if (!(t == u)) {
    ostringstream os;
    os << "Assertion failed: " << t << " != " << u;
    if (!hint.empty()) {
       os << " hint: " << hint;
    }
    throw runtime_error(os.str());
  }
}

template<class T, class U>
void ContainersAssertEqualAccurateTo(const double accuracy, const T& t, const U& u, const string& hint = {}) {
  if (arma::any(arma::vectorise(2 * abs(t - u) / abs(t + u)) > accuracy)) {
    ostringstream os;
    os << "Assertion failed: " << t << " != " << u << ", with accuracy level: " << accuracy;
    if (!hint.empty()) {
       os << " hint: " << hint;
    }
    throw runtime_error(os.str());
  }
}

void DoublesAssertEqualAccurateTo(const double accuracy, const double t, const double u, const string& hint = {}) {
  if ( 2 * abs(t - u) / abs(t + u) > accuracy ) {
    ostringstream os;
    os << "Assertion failed: " << t << " != " << u << ", with accuracy level: " << accuracy;
    if (!hint.empty()) {
       os << " hint: " << hint;
    }
    throw runtime_error(os.str());
  }
}

inline void Assert(bool b, const string& hint) {
  AssertEqual(b, true, hint);
}

class TestRunner {
public:
  template <class TestFunc>
  void RunTest(TestFunc func, const string& test_name) {
    try {
      func();
      cerr << test_name << " OK" << endl;
    } catch (exception& e) {
      ++fail_count;
      cerr << test_name << " fail: " << e.what() << endl;
    } catch (...) {
      ++fail_count;
      cerr << "Unknown exception caught" << endl;
    }
  }

  ~TestRunner() {
    if (fail_count > 0) {
      cerr << fail_count << " unit tests failed. Terminate" << endl;
      exit(1);
    }
  }

private:
  int fail_count = 0;
};

#define DOUBLES_ASSERT_EQUAL_ACCURATE_TO(accuracy, x, y) {               \
  ostringstream os;                                                      \
  os << #x << " != " << #y << ", with accuracy: " <<                     \
    #accuracy << ", " << __FILE__ << ":" << __LINE__;                    \
    DoublesAssertEqualAccurateTo(accuracy, x, y, os.str());              \
}

#define CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(accuracy, x, y) {            \
  ostringstream os;                                                      \
  os << #x << " != " << #y << ", with accuracy: " <<                     \
    #accuracy << ", " << __FILE__ << ":" << __LINE__;                    \
    ContainersAssertEqualAccurateTo(accuracy, x, y, os.str());           \
}

#define ASSERT_EQUAL(x, y) {            \
  ostringstream os;                     \
  os << #x << " != " << #y << ", "      \
    << __FILE__ << ":" << __LINE__;     \
  AssertEqual(x, y, os.str());          \
}

#define ASSERT(x) {                     \
  ostringstream os;                     \
  os << #x << " is false, "             \
    << __FILE__ << ":" << __LINE__;     \
  Assert(x, os.str());                  \
}

#define RUN_TEST(tr, func) \
  tr.RunTest(func, #func)

