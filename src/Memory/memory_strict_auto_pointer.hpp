/// strict_auto_ptr class
template<class T>
class strict_auto_ptr : public std::auto_ptr<T> {
public:
  strict_auto_ptr(T* p = NULL) throw() : std::auto_ptr<T>(p) { }
private:
  strict_auto_ptr (const strict_auto_ptr&) throw();
  void operator = ( const strict_auto_ptr&) throw();
};
