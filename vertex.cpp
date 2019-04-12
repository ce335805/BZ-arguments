#include "vertex.h"

Vertex::Vertex(const int Nl)
    : Nl_internal_ { Nl }
    , vertex_ { std::vector<std::complex<double>>((Nl * 2 * Nw) * (Nl * 2 * Nw) * (NxIBZ * Nw), std::complex<double>(0.0, 0.0)) }
{
}

Vertex::Vertex(std::complex<double> initial_value, const int Nl)
    : Nl_internal_ { Nl }
    , vertex_ { std::vector<std::complex<double>>((Nl * 2 * Nw) * (Nl * 2 * Nw) * (NxIBZ * Nw), initial_value) }
{
}
