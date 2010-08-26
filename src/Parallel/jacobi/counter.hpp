#ifndef COUNTER_HPP
#define COUNTER_HPP

class Counter {
private:
  int count_;
  int count_max_;
public:
  Counter (int count_max)
    : count_max_(count_max),
      count_(count_max)
  {}
  int next()
  {
    count_--;
    if (count_ == 0) count_ = count_max_;
    return count_ == count_max_;
  }
};

#endif
