using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SplineBuilder
{
    public class Spline
    {
        private double[] vertexes;
        private double[] values;
        private double[] coeffs;
        private int size;

        private int i;
        private double h;
        private double ksi;
        private double ksi_2;
        private double ksi_3;

        private double[,] diagonals;
        private static int[] offsets = new int[] { -1, 0, 1 };

        public Spline(double[] vertexes, double[] values)
        {
            if (vertexes.Length != values.Length)
                throw new IndexOutOfRangeException();
            else
            {
                this.vertexes = vertexes;
                this.values = values;
                size = vertexes.Length;
                coeffs = new double[vertexes.Length];
            }

            diagonals = new double[3, vertexes.Length];
        }
        public Spline(List<float[]> points)
        {
            vertexes = new double[points.Count];
            values = new double[points.Count];
            size = vertexes.Length;

            for(int i = 0; i < points.Count; i++)
            {
                vertexes[i] = points[i][0];
                values[i] = points[i][1];
            }
            coeffs = new double[vertexes.Length];
            diagonals = new double[3, vertexes.Length];
        }
        public Spline(double[] vertexes, Func<double, double> values)
        {
            this.vertexes = vertexes;
            size = vertexes.Length;
            coeffs = new double[size];

            this.values = new double[size];
            for (int i = 0; i < size; i++)
                this.values[i] = values(vertexes[i]);

            diagonals = new double[3, size];
        }

        public void CreateSpline()
        {
            double l, r, c;

            double next_h = vertexes[2] - vertexes[1]; // h1
            double prev_h = vertexes[1] - vertexes[0]; // h2

            // f'0
            l = -1.0 * (3.0 * prev_h + 2.0 * next_h) /
                (prev_h * (prev_h + next_h));
            c = (2.0 * next_h + prev_h) /
                (prev_h * next_h);
            r = -1.0 * prev_h /
                (next_h * (prev_h + next_h));
            coeffs[0] = 0.5 * (values[0] * l + values[1] * c + values[2] * r);

            // f'n
            next_h = vertexes[size - 1] - vertexes[size - 2];
            prev_h = vertexes[size - 2] - vertexes[size - 3];
            l = next_h / (prev_h * (prev_h + next_h));
            c = -1.0 * (next_h + 2.0 * prev_h) / (prev_h * next_h);
            r = (3.0 * next_h + 2.0 * prev_h) / (next_h * (prev_h + next_h));
            coeffs[size - 1] = 0.5 * (values[size - 3] * l + values[size - 2] * c + values[size - 1] * r);

            diagonals[1, 0] = 1;
            diagonals[1, size - 1] = 1;

            for (int i = 1; i < size - 1; i++)
            {
                next_h = vertexes[i + 1] - vertexes[i];
                prev_h = vertexes[i] - vertexes[i - 1];

                diagonals[0, i] = 2.0 / prev_h;
                diagonals[1, i] = 4.0 * (1.0 / prev_h + 1.0 / next_h);
                diagonals[2, i] = 2.0 / next_h;

                coeffs[i] = -6.0 * values[i - 1] / (prev_h * prev_h) +
                    6.0 * values[i] * (1.0 / (prev_h * prev_h) - 1.0 / (next_h * next_h)) +
                    6.0 * values[i + 1] / (next_h * next_h);
            }

            solveSeidel(1E-12);
        }

        public int solveSeidel(double epsilon)
        {
            int diagCount = offsets.Length;
            int mainDiagIndex = Array.IndexOf(offsets, 0);
            double[] v = coeffs;
            double[] x = new double[size];
            int i, j;
            double collector;
            int step;
            double norm = 0;
            double curNorm;
            int offset;
            double residual;
            for (i = 0; i < size; i++)
                norm += v[i] * v[i];

            for (step = 0; step < 1000; step++)
            {
                curNorm = 0;
                for (i = 0; i < size; i++)
                {
                    collector = 0;
                    for (j = 0; j < diagCount; j++)
                    {
                        offset = i + offsets[j];
                        if (offset < 0)
                            continue;
                        if (offset >= size)
                            break;
                        collector += diagonals[j, i] * x[offset];
                    }
                    x[i] += 1.5 * (v[i] - collector) / diagonals[mainDiagIndex, i];
                    curNorm += (v[i] - collector) * (v[i] - collector);
                }
                residual = Math.Sqrt(curNorm / norm);
                if (residual < epsilon)
                    break;
            }
            coeffs = x;
            return step + 1;
        }
        public double getValue(double x)
        {
            for (i = 0; i < size - 1; i++)
                if (vertexes[i] <= x && x <= vertexes[i + 1])
                    break;
            if (i == size - 1)
                return double.NaN;

            h = vertexes[i + 1] - vertexes[i];
            ksi = (x - vertexes[i]) / h;
            ksi_2 = ksi * ksi;
            ksi_3 = ksi_2 * ksi;

            return values[i] * (1.0 - 3.0 * ksi_2 + 2.0 * ksi_3) +
                coeffs[i] * h * (ksi - 2.0 * ksi_2 + ksi_3) +
                values[i + 1] * (3.0 * ksi_2 - 2.0 * ksi_3) +
                coeffs[i + 1] * h * (-1.0 * ksi_2 + ksi_3);
        }
    }
}
