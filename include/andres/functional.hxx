#pragma once
#ifndef ANDRES_FUNCTIONAL_HXX
#define ANDRES_FUNCTIONAL_HXX

namespace andres {

template<typename T>
struct Identity {
    typedef T argument_type;
    typedef T result_type;
    constexpr T operator()(const T& x) const
        { return x; }
};

}

#endif
