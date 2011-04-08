#include "test.decl.h"
#include "hello.hpp"
#include "main.hpp"

CProxy_Main  main_proxy;
CProxy_Hello hello_proxy;

Main::Main(CkArgMsg* main) : count_(0)
{
  hello_proxy = CProxy_Hello::ckNew();
}

void Main::p_print ()
{
  if (++count_ == CkNumPes()) {
    Hello * hello = hello_proxy.ckLocal();
    hello->p_print();
    count_ = 0;
  }
}

void Main::p_exit()
{
  CkExit();
}

#include "test.def.h"
