/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        // Calculate weights once instead of 2 times per pixel
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

        const int width = m.get_x_size();
        const int height = m.get_y_size();

        // Save memory dont always create at max size
        Matrix scratch{static_cast<unsigned>(width), static_cast<unsigned>(height)};

        const int leftEnd = (radius < width) ? radius : width;
        const int rightStart = (radius > width - radius) ? radius : width - radius;

        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                double red = m.r(x, y) * w[0];
                double green = m.g(x, y) * w[0];
                double blue = m.b(x, y) * w[0];
                double n = w[0];

                if (x < leftEnd || x >= rightStart)
                {
                    // pixels near edge require bounds check
                    for (int wi = 1; wi <= radius; wi++)
                    {
                        if (x - wi >= 0)
                        {
                            red += m.r(x - wi, y) * w[wi];
                            green += m.g(x - wi, y) * w[wi];
                            blue += m.b(x - wi, y) * w[wi];
                            n += w[wi];
                        }

                        if (x + wi < width)
                        {
                            red += m.r(x + wi, y) * w[wi];
                            green += m.g(x + wi, y) * w[wi];
                            blue += m.b(x + wi, y) * w[wi];
                            n += w[wi];
                        }
                    }
                }
                else
                {
                    // pixels in middle can skip bounds check
                    for (int wi = 1; wi <= radius; wi++)
                    {
                        red += m.r(x - wi, y) * w[wi];
                        green += m.g(x - wi, y) * w[wi];
                        blue += m.b(x - wi, y) * w[wi];
                        n += w[wi];

                        red += m.r(x + wi, y) * w[wi];
                        green += m.g(x + wi, y) * w[wi];
                        blue += m.b(x + wi, y) * w[wi];
                        n += w[wi];
                    }
                }

                scratch.r(x, y) = red / n;
                scratch.g(x, y) = green / n;
                scratch.b(x, y) = blue / n;
            }
        }

        const int topEnd = (radius < height) ? radius : height;
        const int bottomStart = (radius > height - radius) ? radius : height - radius;

        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                double red = scratch.r(x, y) * w[0];
                double green = scratch.g(x, y) * w[0];
                double blue = scratch.b(x, y) * w[0];
                double n = w[0];

                if (y < topEnd || y >= bottomStart)
                {
                    // pixels near edge require bounds check
                    for (int wi = 1; wi <= radius; wi++)
                    {
                        if (y - wi >= 0)
                        {
                            red += scratch.r(x, y - wi) * w[wi];
                            green += scratch.g(x, y - wi) * w[wi];
                            blue += scratch.b(x, y - wi) * w[wi];
                            n += w[wi];
                        }

                        if (y + wi < height)
                        {
                            red += scratch.r(x, y + wi) * w[wi];
                            green += scratch.g(x, y + wi) * w[wi];
                            blue += scratch.b(x, y + wi) * w[wi];
                            n += w[wi];
                        }
                    }
                }
                else
                {
                    // pixels in middle can skip bounds check
                    for (int wi = 1; wi <= radius; wi++)
                    {
                        red += scratch.r(x, y - wi) * w[wi];
                        green += scratch.g(x, y - wi) * w[wi];
                        blue += scratch.b(x, y - wi) * w[wi];
                        n += w[wi];

                        red += scratch.r(x, y + wi) * w[wi];
                        green += scratch.g(x, y + wi) * w[wi];
                        blue += scratch.b(x, y + wi) * w[wi];
                        n += w[wi];
                    }
                }

                m.r(x, y) = red / n;
                m.g(x, y) = green / n;
                m.b(x, y) = blue / n;
            }
        }

        return m;
    }
}
