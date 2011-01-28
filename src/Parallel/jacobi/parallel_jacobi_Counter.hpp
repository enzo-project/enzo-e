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
	count_(0)
    {}
    int remaining()
    {
      return count_ = (count_max_ + count_ - 1) % count_max_;
    }
  };
}
#endif /* PARALLEL_JACOBI_COUNTER_HPP */

