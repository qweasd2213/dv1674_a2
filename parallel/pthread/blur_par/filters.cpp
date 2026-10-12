/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <vector>
#include <pthread.h>
#include <algorithm>

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

    struct BlurTask
    {
        Matrix *src;
        Matrix *dst;
        const double *weights;
        int radius;
        int startX;
        int endX;
        bool horizontal;
    };

    void *blur_thread(void *arg)
    {
        auto *task = static_cast<BlurTask *>(arg);

        Matrix &src = *task->src;
        Matrix &dst = *task->dst;
        const double *w = task->weights;

        const int width = src.get_x_size();
        const int height = src.get_y_size();
        const int radius = task->radius;

        if (task->horizontal)
        {
            const int leftEnd = (radius < width) ? radius : width;
            const int rightStart = (radius > width - radius) ? radius : width - radius;

            for (int x = task->startX; x < task->endX; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    double red = src.r(x, y) * w[0];
                    double green = src.g(x, y) * w[0];
                    double blue = src.b(x, y) * w[0];
                    double n = w[0];

                    if (x < leftEnd || x >= rightStart)
                    {
                        for (int wi = 1; wi <= radius; wi++)
                        {
                            if (x - wi >= 0)
                            {
                                red += src.r(x - wi, y) * w[wi];
                                green += src.g(x - wi, y) * w[wi];
                                blue += src.b(x - wi, y) * w[wi];
                                n += w[wi];
                            }

                            if (x + wi < width)
                            {
                                red += src.r(x + wi, y) * w[wi];
                                green += src.g(x + wi, y) * w[wi];
                                blue += src.b(x + wi, y) * w[wi];
                                n += w[wi];
                            }
                        }
                    }
                    else
                    {
                        for (int wi = 1; wi <= radius; wi++)
                        {
                            red += src.r(x - wi, y) * w[wi];
                            green += src.g(x - wi, y) * w[wi];
                            blue += src.b(x - wi, y) * w[wi];
                            n += w[wi];

                            red += src.r(x + wi, y) * w[wi];
                            green += src.g(x + wi, y) * w[wi];
                            blue += src.b(x + wi, y) * w[wi];
                            n += w[wi];
                        }
                    }

                    dst.r(x, y) = red / n;
                    dst.g(x, y) = green / n;
                    dst.b(x, y) = blue / n;
                }
            }
        }
        else
        {
            const int topEnd = (radius < height) ? radius : height;
            const int bottomStart =
                (radius > height - radius) ? radius : height - radius;

            for (int x = task->startX; x < task->endX; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    double red = src.r(x, y) * w[0];
                    double green = src.g(x, y) * w[0];
                    double blue = src.b(x, y) * w[0];
                    double n = w[0];

                    if (y < topEnd || y >= bottomStart)
                    {
                        for (int wi = 1; wi <= radius; wi++)
                        {
                            if (y - wi >= 0)
                            {
                                red += src.r(x, y - wi) * w[wi];
                                green += src.g(x, y - wi) * w[wi];
                                blue += src.b(x, y - wi) * w[wi];
                                n += w[wi];
                            }

                            if (y + wi < height)
                            {
                                red += src.r(x, y + wi) * w[wi];
                                green += src.g(x, y + wi) * w[wi];
                                blue += src.b(x, y + wi) * w[wi];
                                n += w[wi];
                            }
                        }
                    }
                    else
                    {
                        for (int wi = 1; wi <= radius; wi++)
                        {
                            red += src.r(x, y - wi) * w[wi];
                            green += src.g(x, y - wi) * w[wi];
                            blue += src.b(x, y - wi) * w[wi];
                            n += w[wi];

                            red += src.r(x, y + wi) * w[wi];
                            green += src.g(x, y + wi) * w[wi];
                            blue += src.b(x, y + wi) * w[wi];
                            n += w[wi];
                        }
                    }

                    dst.r(x, y) = red / n;
                    dst.g(x, y) = green / n;
                    dst.b(x, y) = blue / n;
                }
            }
        }

        return nullptr;
    }

    Matrix blur(Matrix m, const int radius, int threads)
    {        
        // Calculate weights once instead of 2 times per pixel
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

        const int width = m.get_x_size();
        const int height = m.get_y_size();

        // Save memory dont always create at max size
        Matrix scratch{static_cast<unsigned>(width), static_cast<unsigned>(height)};

        std::vector<pthread_t> threadIds(threads);
        std::vector<BlurTask> tasks(threads);

        const int columnsPerThread = width / threads;
        const int remainder = width % threads;

        int startX = 0;

        for (int i = 0; i < threads; i++)
        {
            const int count = columnsPerThread + (i < remainder ? 1 : 0);

            tasks[i] = {
                &m, &scratch, w, radius,
                startX, startX + count, true
            };

            pthread_create(
                &threadIds[i], nullptr, blur_thread, &tasks[i]);

            startX += count;
        }

        for (int i = 0; i < threads; i++)
            pthread_join(threadIds[i], nullptr);

        startX = 0;

        for (int i = 0; i < threads; i++)
        {
            const int count = columnsPerThread + (i < remainder ? 1 : 0);

            tasks[i] = {
                &scratch, &m, w, radius,
                startX, startX + count, false
            };

            pthread_create(
                &threadIds[i], nullptr, blur_thread, &tasks[i]);

            startX += count;
        }

        for (int i = 0; i < threads; i++)
            pthread_join(threadIds[i], nullptr);

        return m;
    }
}
