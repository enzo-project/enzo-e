

//----------------------------------------------------------------------
#include "test.decl.h"
//----------------------------------------------------------------------

class GroupBASE { 
public:
  GroupBASE () { CkPrintf ("GroupBASE()\n");  };
  virtual ~GroupBASE() {}; 
};

class ArrayBASE { 
public:
  ArrayBASE (int n) { CkPrintf ("ArrayBASE(%d)\n",n); };
  ArrayBASE (CkMigrateMessage *) {};
  ArrayBASE () {};
  virtual ~ArrayBASE() {}; 
};


class HelloGroup : public Group, public GroupBASE {
public:
  HelloGroup() ;
  ~HelloGroup() {delete [] print_; };
  void p_print();
private:
  char * print_;
   
};

//----------------------------------------------------------------------

class HelloArray : public CBase_HelloArray , public ArrayBASE {
public:
  HelloArray() {};
  HelloArray(int n) ;
  HelloArray (CkMigrateMessage *) {};
  ~HelloArray() {delete [] print_; };
  void p_print();
private:
  char * print_;
};

//----------------------------------------------------------------------

class Main : public Chare {

public:
	   
  Main(CkArgMsg *m);  // Calls HelloGroup()
  void create_group();
  void p_print_group ();  // Calls HelloGroup p_print()
  void p_create_array();  // Creates new HelloArray
  void p_print_array ();  // Calls HelloArray p_print()
  void p_exit(); // Exit at end

private:
  int count_;
    
};

//======================================================================

CProxy_Main       main_proxy;
CProxy_HelloGroup hello_group_proxy;
CProxy_HelloArray hello_array_proxy;

//----------------------------------------------------------------------

Main::Main(CkArgMsg* main) : count_(0)
{
  create_group();
}

//----------------------------------------------------------------------

void Main::create_group ()
{
  // Create group, switching from Main to Group
  hello_group_proxy = CProxy_HelloGroup::ckNew();
}

//----------------------------------------------------------------------

HelloGroup::HelloGroup() 
{ 
  print_ = new char[80];
  sprintf (print_, "HelloGroup node %d proc %d\n",CkMyNode(),CkMyPe());
  main_proxy.p_print_group();
};

//----------------------------------------------------------------------

void Main::p_print_group ()
{
  // called by HelloGroup
  if (++count_ == CkNumPes()) {
    // Main
    hello_group_proxy.p_print();
    count_ = 0;
  }
}

//----------------------------------------------------------------------

void HelloGroup::p_print()
  {
    CkPrintf ("%s",print_);
    main_proxy.p_create_array();
  };


//----------------------------------------------------------------------

void Main::p_create_array()
{
  // Called by HelloGroup::p_print()
  if (++count_ == CkNumPes()) {
    hello_array_proxy = CProxy_HelloArray::ckNew(10,10);
  }
}

//----------------------------------------------------------------------

HelloArray::HelloArray(int n) 
{ 
  print_ = new char[80];
  sprintf (print_, "HelloArray element %d/%d\n",thisIndex,n);
  main_proxy.p_print_array();
};

//----------------------------------------------------------------------

void Main::p_print_array ()
{
  if (++count_ == 10) {
    hello_array_proxy.p_print();
    count_ = 0;
  }
}

//----------------------------------------------------------------------

void HelloArray::p_print()
{
  CkPrintf ("%s",print_);
  main_proxy.p_exit();
};

//----------------------------------------------------------------------

void Main::p_exit()
{
  if (++count_ == 10) {
    CkPrintf ("Exiting\n");
    CkExit();
  }
}

//----------------------------------------------------------------------

#include "test.def.h"
