#ifndef CommonTools_Utils_EtSCComparator_h
#define CommonTools_Utils_EtSCComparator_h

template<typename T>
struct LessBySCEt {
  typedef T first_argument_type;
  typedef T second_argument_type;
  bool operator()( const T & t1, const T & t2 ) const {
    float et1 = t1.energy() * sin( t1.position().theta() );
    float et2 = t2.energy() * sin( t2.position().theta() );
    return et1 < et2;
  }
};

template<typename T>
struct GreaterBySCEt {
  typedef T first_argument_type;
  typedef T second_argument_type;
  bool operator()( const T & t1, const T & t2 ) const {
    float et1 = t1.energy() * sin( t1.position().theta() );
    float et2 = t2.energy() * sin( t2.position().theta() );
    return et1 > et2;
  }
};


#endif
