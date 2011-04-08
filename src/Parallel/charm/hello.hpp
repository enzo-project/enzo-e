
extern CProxy_Main main_proxy;

class Hello : public Group {
public:
  Hello() 
  { 
    print_ = new char[80];
    sprintf (print_, "Hello node %d proc %d\n",CkMyNode(),CkMyPe());
    main_proxy.p_print();
  };

  void p_print()
  {
    printf ("%s",print_);
    main_proxy.p_exit();
  };

private:
  char * print_;
    
};
