#pragma once
#ifndef ANDRES_GRAPH_DO_NOTING_FUNCTOR_HXX
#define ANDRES_GRAPH_DO_NOTING_FUNCTOR_HXX

namespace andres {
namespace graph {
namespace detail {

template<typename T>
struct do_nothing
{
    typedef T argument_type;
    typedef T result_type;

    constexpr T operator()(const T& x) const
    {
        return x;
    }
};

}
}
}
#endif
