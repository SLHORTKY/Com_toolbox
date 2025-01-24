#pragma once
#include "Array.h"

namespace Com
{
    class Sig
    {
    private:
        Array x; 
        Array y;
        
    public:
        Sig(Array X, Array Y);
        ~Sig() = default;

        Array getX();
        Array getY();
    };
    
    inline Array Sig::getX(){return this->x;}
    inline Array Sig::getY(){return this->y;}

    inline Sig::Sig(Array X, Array Y): x(X),y(Y){}
} // namespace Com
