#pragma once
#ifndef ANDRES_RUNTIME_CHECK_HXX
#define ANDRES_RUNTIME_CHECK_HXX

#include <cstdint>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>




/** \def GRAPH_CHECK_OP(a,op,b,message)
    \brief macro for runtime checks
    
    \warning The check is done 
        <B> even in Release mode </B> 
        (therefore if NDEBUG <B>is</B> defined)

    \param a : first argument (like a number )
    \param op : operator (== )
    \param b : second argument (like a number )
    \param message : error message (as "my error")

    <b>Usage:</b>
    \code
        int a = 1;
        GRAPH_CHECK_OP(a, ==, 1, "this should never fail")
        GRAPH_CHECK_OP(a, >=, 2, "this should fail")
    \endcode
*/
#define GRAPH_CHECK_OP(a,op,b,message) \
    if(!  static_cast<bool>( a op b )   ) { \
       std::stringstream s; \
       s << "Graph Error: "<< message <<"\n";\
       s << "Graph check :  " << #a <<#op <<#b<< "  failed:\n"; \
       s << #a " = "<<a<<"\n"; \
       s << #b " = "<<b<<"\n"; \
       s << "in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
       throw std::runtime_error(s.str()); \
    }

/** \def GRAPH_CHECK(expression,message)
    \brief macro for runtime checks
    
    \warning The check is done 
        <B> even in Release mode </B> 
        (therefore if NDEBUG <B>is</B> defined)

    \param expression : expression which can evaluate to bool
    \param message : error message (as "my error")

    <b>Usage:</b>
    \code
        int a = 1;
        GRAPH_CHECK_OP(a==1, "this should never fail")
        GRAPH_CHECK_OP(a>=2, "this should fail")
    \endcode
*/
#define GRAPH_CHECK(expression, message) if(!(expression)) { \
   std::stringstream s; \
   s << message <<"\n";\
   s << "Graph assertion " << #expression \
   << " failed in file " << __FILE__ \
   << ", line " << __LINE__ << std::endl; \
   throw std::runtime_error(s.str()); \
 }


#define GRAPH_TEST(expression) GRAPH_CHECK(expression,"")

#define GRAPH_TEST_OP(a,op,b) GRAPH_CHECK_OP(a,op,b,"")

#define GRAPH_CHECK_EQ_TOL(a,b,tol) \
    if( std::abs(a-b) > tol){ \
        std::stringstream s; \
        s<<"Graph assertion "; \
        s<<"\""; \
        s<<" | "<< #a <<" - "<<#b<<"| < " #tol <<"\" "; \
        s<<"  failed with:\n"; \
        s<<#a<<" = "<<a<<"\n";\
        s<<#b<<" = "<<b<<"\n";\
        s<<#tol<<" = "<<tol<<"\n";\
        throw std::runtime_error(s.str()); \
    }

#define GRAPH_CHECK_NUMBER(number) \
   { \
   std::stringstream s; \
   s << "Graph assertion failed in file " << __FILE__ \
   << ", line " << __LINE__ << std::endl; \
    if(std::isnan(number))\
        throw std::runtime_error(s.str()+" number is nan"); \
    if(std::isinf(number))\
        throw std::runtime_error(s.str()+"number is inf");\
    }

#ifdef NDEBUG
    #ifdef GRAPH_DEBUG 
        #define GRAPH_DO_DEBUG
    #endif
#else
    #ifdef GRAPH_DEBUG 
        #define GRAPH_DO_DEBUG
    #endif
#endif


/** \def GRAPH_ASSERT_OP(a,op,b,message)
    \brief macro for runtime checks
    
    \warning The check is <B>only</B> done in
        in Debug mode (therefore if NDEBUG is <B>not</B> defined)

    \param a : first argument (like a number )
    \param op : operator (== )
    \param b : second argument (like a number )
    \param message : error message (as "my error")

    <b>Usage:</b>
    \code
        int a = 1;
        GRAPH_ASSERT_OP(a, ==, 1) // will not fail here
        GRAPH_ASSERT_OP(a, >=, 2) // will fail here
    \endcode
*/
#ifdef NDEBUG
   #ifndef GRAPH_DEBUG 
      #define GRAPH_ASSERT_OP(a,op,b) { }
   #else
      #define GRAPH_ASSERT_OP(a,op,b) \
      if(!  static_cast<bool>( a op b )   ) { \
         std::stringstream s; \
         s << "Graph assertion :  " << #a <<#op <<#b<< "  failed:\n"; \
         s << #a " = "<<a<<"\n"; \
         s << #b " = "<<b<<"\n"; \
         s << "in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
         throw std::runtime_error(s.str()); \
      }
   #endif
#else
   #define GRAPH_ASSERT_OP(a,op,b) \
   if(!  static_cast<bool>( a op b )   ) { \
      std::stringstream s; \
      s << "Graph assertion :  " << #a <<#op <<#b<< "  failed:\n"; \
      s << #a " = "<<a<<"\n"; \
      s << #b " = "<<b<<"\n"; \
      s << "in file " << __FILE__ << ", line " << __LINE__ << "\n"; \
      throw std::runtime_error(s.str()); \
   }
#endif

/** \def GRAPH_ASSERT(expression,message)
    \brief macro for runtime checks
    
    \warning The check is <B>only</B> done in
        in Debug mode (therefore if NDEBUG is <B>not</B> defined)

    \param expression : expression which can evaluate to bool

    <b>Usage:</b>
    \code
        int a = 1;
        GRAPH_ASSERT(a == 1) // will not fail here 
        GRAPH_ASSERT(a >= 2) // will fail here
    \endcode
*/
#ifdef NDEBUG
   #ifndef GRAPH_DEBUG
      #define GRAPH_ASSERT(expression) {}
   #else
      #define GRAPH_ASSERT(expression) if(!(expression)) { \
         std::stringstream s; \
         s << "Graph assertion " << #expression \
         << " failed in file " << __FILE__ \
         << ", line " << __LINE__ << std::endl; \
         throw std::runtime_error(s.str()); \
      }
   #endif
#else
      #define GRAPH_ASSERT(expression) if(!(expression)) { \
         std::stringstream s; \
         s << "Graph assertion " << #expression \
         << " failed in file " << __FILE__ \
         << ", line " << __LINE__ << std::endl; \
         throw std::runtime_error(s.str()); \
      }
#endif




#endif // #ifndef ANDRES_RUNTIME_CHECK_HXX
