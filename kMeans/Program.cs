using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Clustering;
using System.Threading.Tasks;

namespace kMeans
{
     class Program
    {
        static void Main(string[] args)
        {
            Random r = new Random();
            List<double> X = new List<double>(400); List<double> Y = new List<double>(400);
            double centroid1X = 4; double centroid1Y = 4;
            double centroid2X = -4; double centroid2Y = -4;
            double centroid3X = -4; double centroid3Y = 4;
            double centroid4X = 4; double centroid4Y = -4;

            double mean = 0;
            double stdDev = 1.4;
            for(int i = 0; i<100; i++)
            {
                X.Add(centroid1X + GaussianNoise(r, mean, stdDev));
                X.Add(centroid2X + GaussianNoise(r, mean, stdDev));
                X.Add(centroid3X + GaussianNoise(r, mean, stdDev));
                X.Add(centroid4X + GaussianNoise(r, mean, stdDev));
                Y.Add(centroid1Y + GaussianNoise(r, mean, stdDev));
                Y.Add(centroid2Y + GaussianNoise(r, mean, stdDev));
                Y.Add(centroid3Y + GaussianNoise(r, mean, stdDev));
                Y.Add(centroid4Y + GaussianNoise(r, mean, stdDev));
            }

            List<Coordinates> Centroids = Clustering.kMeans.GetCentroids(X, Y, 4, 50);

            for(int i = 0; i<4; i++)
            {
                Console.WriteLine(Centroids.ElementAt(i).x.ToString() + " " + Centroids.ElementAt(i).y.ToString());
                
            }
            Console.Read();
        }

        public static double GaussianNoise(Random rand, double mean, double stdDev)
        {
            double u1 = 1.0 - rand.NextDouble(); //uniform(0,1] random doubles
            double u2 = 1.0 - rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            double randNormal =
                         mean + stdDev * randStdNormal; //random normal(mean,stdDev^2)
            return randNormal;
        }
    }
}
