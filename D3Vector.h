#ifndef _D3Vector_
#define _D3Vector_ 

template< class T >
class D3Vector {
 
template< class U >
friend std::ostream & operator<<( std::ostream & , const D3Vector<U> & ) ;   
 
public :
   D3Vector( T a , T b , T c ) {
      x = a ;
      y = b ;
      z = c ;
      w=1.0f;
    }

  D3Vector addi ( const D3Vector & rhs ) {
    
      T a = x + rhs.x;  
      T b = y + rhs.y;  
      T c = z + rhs.z;  
    D3Vector add(a,b,c);
    return add;
   }
  
  D3Vector subs ( const D3Vector & rhs ) {
    
      T a = x - rhs.x;  
      T b = y - rhs.y;  
      T c = z - rhs.z;  
    D3Vector sub(a,b,c);
    return sub;
   }
  
  D3Vector scalarp ( double rhs ) {
    
      T a = x * rhs;  
      T b = y * rhs;  
      T c = z * rhs;  
    D3Vector scalp(a,b,c);
    return scalp;
   }

   T dotp ( const D3Vector & rhs ) {
      T scalar = x * rhs.x + y * rhs.y + z * rhs.z ;
      return scalar ;
   }
 
   D3Vector crossp ( const D3Vector & rhs ) {
      T a = y * rhs.z - z * rhs.y ;
      T b = z * rhs.x - x * rhs.z ;
      T c = x * rhs.y - y * rhs.x ;
      D3Vector product( a , b , c ) ;
      return product ;
   }

  double copyX() {
   T xcomp = x;
  return xcomp;
   }
  double copyY() {
   T ycomp = y;
  return ycomp;
   }
  double copyZ() {
   T zcomp = z;
  return zcomp;
   }

   double length () {
      T len = sqrt(x * x + y * y + z * z);
      return len ;
   }
   
   void normalize()
   {
        double len = length();
        x = x/len;
        y = y/len;
        z = z/len;
        w =1.0;
   }

   private :
   T x , y , z , w;  
 };
 
template< class T >
std::ostream & operator<< ( std::ostream & os , const D3Vector<T> & vec ) {
   os << "( "  << vec.x << " ,  " << vec.y << " ,  " << vec.z << " )" ;
   return os ;
}
 
#endif
