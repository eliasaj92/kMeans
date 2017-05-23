using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Clustering
{
    public class Coordinates
    {
        public double x;
        public double y;
    }

    class kMeans
    {
        public static List<Coordinates> GetCentroids(List<double> X, List<double> Y, int numOfCentroids, int numOfIterations = 20, double tolerableDiff = 0)
        {
            List<Coordinates> Centroids = new List<Coordinates>();
            double Xmean, Xfactor, Ymean, Yfactor;
            Normalize(ref X, out Xmean, out Xfactor);
            Normalize(ref Y, out Ymean, out Yfactor);
            Centroids = InitializeCentroids(X, Y, numOfCentroids);
            List<int> inCluster = new List<int>();
            double diffBetweenIterations = double.PositiveInfinity;
            for (int i = 0; i < numOfIterations; i++)
            {
                inCluster.Clear();
                UpdateClusters(ref inCluster, X, Y, Centroids);

                UpdateCentroids(ref Centroids, X, Y, numOfCentroids, inCluster, out diffBetweenIterations);
                if(diffBetweenIterations<tolerableDiff)
                {
                    break;
                }
            }

            Denormalize(ref Centroids, Xmean, Xfactor, Ymean, Yfactor);
            return Centroids;
        }

        public static List<Coordinates> InitializeCentroids(List<double> X, List<double> Y, int numOfCentroids)
        {
            /*instead of a random initialization like the usual k-means algorithm, this one puts in the initial cluster centers in a quadrilateral whose edges
            are based on the mins and maxes of either coordinate. This is a specific implementation for better clustering for recatangular QAM symbol maps
            */
            double minX = X.Min(); double maxX = X.Max();
            double minY = Y.Min(); double maxY = Y.Max();
            double diffX = maxX - minX; double diffY = maxY - minY;

            List<Coordinates> InitialPoints = new List<Coordinates>();

            for (int i = 0; i < Math.Sqrt(numOfCentroids); i++)
            {
                double tempY = minY + diffY / (Math.Sqrt(numOfCentroids) - 1) * (i % (Math.Round(Math.Sqrt(numOfCentroids))));
                for (int j = 0; j < Math.Sqrt(numOfCentroids); j++)
                {
                    double tempX = minX + diffX / (Math.Sqrt(numOfCentroids) - 1) * (j % (Math.Round(Math.Sqrt(numOfCentroids))));
                    Coordinates temp = new Coordinates();
                    temp.x = tempX;
                    temp.y = tempY;
                    InitialPoints.Add(temp);
                }
            }
            //from here on the code does some tweeking in case the nuber of points is not a square integer. Otherwise, the initial points have been placed to form the quadrilateral 

            if (InitialPoints.Count < numOfCentroids)
            {
                int cnt = InitialPoints.Count;
                for (int i = 0; i < numOfCentroids - cnt; i++)
                {
                    Random rnd = new Random();
                    double tempX = rnd.Next((int)minX, (int)maxX);
                    double tempY = rnd.Next((int)minY, (int)maxY);
                    Coordinates temp = new Coordinates();
                    temp.x = tempX;
                    temp.y = tempY;
                    InitialPoints.Add(temp);
                }
            }
            else if (InitialPoints.Count > numOfCentroids)
            {
                int cnt = InitialPoints.Count;
                for (int i = 0; i < cnt - numOfCentroids; i++)
                {
                    InitialPoints.RemoveAt(InitialPoints.Count - 1);
                }
            }

            return InitialPoints;

        }

        public static void UpdateClusters(ref List<int> inCluster, List<double> X, List<double> Y, List<Coordinates> Centroids)
        {
            //Classifies points into clusters based on nearest centroid
            //inCluster is a list whose elements are the cluster number and whose index corresponds to the ones in X and Y
            //for example inCluster.AtElement(0) = 1 means that the point with X.AtElement(0) and Y.AtElement(0) belongs to the Cluster 1
            inCluster.Clear();
            int cnt = Math.Min(X.Count, Y.Count);
            for (int i = 0; i < cnt; i++)
            {
                int minIndex = -1;
                double minDist = double.PositiveInfinity;
                for (int j = 0; j < Centroids.Count; j++)
                {
                    double dist = (Centroids.ElementAt(j).x - X.ElementAt(i)) * (Centroids.ElementAt(j).x - X.ElementAt(i))
                        + (Centroids.ElementAt(j).y - Y.ElementAt(i)) * (Centroids.ElementAt(j).y - Y.ElementAt(i));
                    if (dist <= minDist)
                    {
                        minDist = dist;
                        minIndex = j;
                    }
                }
                inCluster.Add(minIndex);
            }

        }

        public static void UpdateCentroids(ref List<Coordinates> Centroids, List<double> X, List<double> Y, int numOfCentroids, List<int> inCluster, out double diff)
        {
            //gets new centroids based on the classified points
            List<Coordinates> tempCentroids = new List<Coordinates>();
            double[] accumulatedX = new double[numOfCentroids];
            double[] accumulatedY = new double[numOfCentroids];
            int[] clusterCount = new int[numOfCentroids];

            for (int j = 0; j < numOfCentroids; j++)
            {
                accumulatedX[j] = 0.0;
                accumulatedY[j] = 0.0;
                clusterCount[j] = inCluster.Where(c => c == j).Count();
            }

            int cnt = Math.Min(X.Count, Y.Count);
            for (int i = 0; i < cnt; i++)
            {
                accumulatedX[inCluster.ElementAt(i)] += X.ElementAt(i);
                accumulatedY[inCluster.ElementAt(i)] += Y.ElementAt(i);
            }
            
            for (int j = 0; j < numOfCentroids; j++)
            {
                Coordinates temp = new Coordinates();
                temp.x = accumulatedX[j] / ((double)clusterCount[j]);
                temp.y = accumulatedY[j] / ((double)clusterCount[j]);
                tempCentroids.Add(temp);
            }

            diff = 0;
            for (int j = 0; j < tempCentroids.Count; j++)
            {
                diff += Math.Sqrt((Centroids.ElementAt(j).x - tempCentroids.ElementAt(j).x) * (Centroids.ElementAt(j).x - tempCentroids.ElementAt(j).x)
                    + (Centroids.ElementAt(j).y - tempCentroids.ElementAt(j).y) * (Centroids.ElementAt(j).y - tempCentroids.ElementAt(j).y));

                //the following is in case one of the centroids was too far from one of the clusters, hence when computing the new centroid
                //a divide by zero was done. So we randomize the centroid location.
                if (Double.IsNaN(tempCentroids.ElementAt(j).x) || Double.IsNaN(tempCentroids.ElementAt(j).y))
                {
                    tempCentroids.RemoveAt(j);
                    Coordinates temp = new Coordinates();
                    double minX = X.Min(); double maxX = X.Max();
                    double minY = Y.Min(); double maxY = Y.Max();
                    Random rnd = new Random();
                    temp.x = rnd.Next((int)minX, (int)maxX);
                    temp.y = rnd.Next((int)minY, (int)maxY);
                    tempCentroids.Add(temp);
                }
            }
            Centroids.Clear();
            Centroids = tempCentroids;
        }

        public static void Normalize(ref List<double>X, out double Xmean, out double Xfactor)
        {
            List<double> tempX = new List<double>();
            double accumulatedX =0.0;
            for(int i = 0; i<X.Count; i++)
            {
                accumulatedX += X.ElementAt(i);
            }

            Xmean = accumulatedX / X.Count;

            for (int i = 0; i < X.Count; i++)
            {
                tempX.Add(X.ElementAt(i) - Xmean);
            }

            double minX = tempX.Min();
            double maxX = tempX.Max();

            Xfactor = maxX - minX;
            X.Clear();

            for(int i = 0; i < tempX.Count; i++)
            {
                X.Add(tempX.ElementAt(i) / Xfactor);
            }
        }

        public static void Denormalize (ref List<Coordinates> Centroids, double Xmean, double Xfactor, double Ymean, double Yfactor )
        {
            List<Coordinates> tempCentroids = new List<Coordinates>();
            for (int j = 0; j < Centroids.Count; j++)
            {
                tempCentroids.Add(Centroids.ElementAt(j));
            }
            Centroids.Clear();
            for(int j = 0; j<tempCentroids.Count; j++)
            {
                Coordinates temp = new Coordinates();
                temp.x = tempCentroids.ElementAt(j).x * Xfactor + Xmean;
                temp.y = tempCentroids.ElementAt(j).y * Yfactor + Ymean;
                Centroids.Add(temp);
            }
        }

    }

}
