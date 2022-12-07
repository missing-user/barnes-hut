#ifndef RECTANGLE
#define RECTANGLE
class Rectangle
{
public:
    double x1, y1, x2, y2;
    Rectangle() {}
    Rectangle(double x1In, double y1In, double x2In, double y2In)
    {
        x1 = x1In;
        y1 = y1In;
        x2 = x2In;
        y2 = y2In;
    }
    bool contains(const Particle &P)
    {
        return P.p.x > x1 && P.p.x < x2 && P.p.y > y1 && P.p.y < y2;
    }
    std::vector<Rectangle> subdivide()
    {
        std::vector<Rectangle> subrectangles(4);
        subrectangles[0] = Rectangle(x1, y1, (x1 + x2) / 2, (y1 + y2) / 2);
        subrectangles[1] = Rectangle((x1 + x2) / 2, y1, x2, (y1 + y2) / 2);
        subrectangles[2] = Rectangle(x1, (y1 + y2) / 2, (x1 + x2) / 2, y2);
        subrectangles[3] = Rectangle((x1 + x2) / 2, (y1 + y2) / 2, x2, y2);
        return subrectangles;
    }
    std::string print() const
    {
        std::ostringstream str;
        str << x1 << " " << y1 << " " << x2 << " " << y2 << " ";
        std::string s = str.str();
        return s;
    }
};

#endif