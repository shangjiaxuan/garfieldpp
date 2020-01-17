#ifndef G_POLYGON_H
#define G_POLYGON_H

#include <vector>

namespace Garfield {

namespace Polygon {

/// Determine whether the point (x, y) is located inside of the
/// polygon (xpl, ypl).
void Inside(const std::vector<double>& xpl, const std::vector<double>& ypl,
            const double x, const double y, bool& inside, bool& edge);

/// Determine the (signed) area of a polygon.
double Area(const std::vector<double>& xp, const std::vector<double>& yp);

}

}

#endif
