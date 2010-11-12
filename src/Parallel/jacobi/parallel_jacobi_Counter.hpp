#ifndef PARALLEL_JACOBI_COUNTER_HPP
#define PARALLEL_JACOBI_COUNTER_HPP

namespace jacobi
{
  class Counter {
  private:
    int count_max_;
    int count_;
  public:
    Counter (int count_max)
      : count_max_(count_max),
	count_(count_max)
    {}
    int wait()
    {
      count_--;
      if (count_ == 0) count_ = count_max_;
      return count_ == count_max_;
    }
  };
}
#endif /* PARALLEL_JACOBI_COUNTER_HPP */

