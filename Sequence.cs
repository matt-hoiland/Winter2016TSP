using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    internal class Sequence : IComparable
    {
        private List<City> route;

        private double score;

        public City this[int i]
        {
            get
            {
                return route[i];
            }
            set
            {
                route[i] = value;
            }
        }

        public int Length
        {
            get
            {
                return route.Count;
            }
        }

        public Sequence(List<City> thing)
        {
            route = thing;
            score = double.PositiveInfinity;
        }

        public Sequence(int[] perm, City[] cities)
        {
            route = new List<City>();
            foreach (int index in perm)
            {
                route.Add(cities[index]);
            }
            Rescore();
        }

        public Sequence(Sequence parent)
        {
            route = new List<City>();
            for (int i = 0; i < parent.route.Count; i++)
            {
                route.Add(parent[i]);
            }
            Rescore();
        }

        public Sequence(ProblemAndSolver.TSPSolution best)
        {
            route = new List<City>();
            foreach (City city in best.Route)
            {
                route.Add(city);
            }
            Rescore();
            /*
            var list = best.Route;
            route = new int[list.Count];
            City[] cities = parent.GetCities();
            int j = 0;
            foreach (City c in list)
            {
                for (int i = 0; i < cities.Length; i++)
                {
                    if (cities[i].Equals(c))
                    {
                        route[j] = i;
                    }
                }
                j++;
            }
            rescore();
            */
        }

        public double Score
        {
            get
            {
                if (double.IsPositiveInfinity(score))
                {
                    Rescore();
                }
                return score;
            }
        }

        public void Rescore()
        {
            score = ToRoute().costOfRoute();
        }

        public int IndexOf(City city)
        {
            return route.IndexOf(city);
            /*
            for (int i = 0; i < Length; i++)
            {
                if (route[i] == city) { return i; }
            }
            return -1;
            */
        }

        public ProblemAndSolver.TSPSolution ToRoute()
        {
            ArrayList list = new ArrayList();
            foreach (var i in route)
            {
                list.Add(i);
            }
            return new ProblemAndSolver.TSPSolution(list);
        }

        public int CompareTo(object obj)
        {
            return Score.CompareTo(((Sequence)obj).Score);
        }
    }
}
