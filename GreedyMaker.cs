using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    class GreedyMaker
    {
        internal class DistanceCity : IComparable
        {
            public double Distance { get; set; }
            public City Destination { get; set; }

            public DistanceCity(double dist, City dest)
            {
                Distance = dist;
                Destination = dest;
            }

            public int CompareTo(object obj)
            {
                return Distance.CompareTo(((DistanceCity)obj).Distance);
            }
        }

        private Dictionary<City, List<DistanceCity>> cities;

        public GreedyMaker(City[] arrCities)
        {
            cities = new Dictionary<City, List<DistanceCity>>();
            foreach (City source in arrCities)
            {
                cities[source] = new List<DistanceCity>();
                foreach (City destination in arrCities)
                {
                    if (destination != source)
                    {
                        cities[source].Add(new DistanceCity(source.costToGetTo(destination), destination));
                    }
                }
                cities[source].Sort();
            }
        }

        internal List<DistanceCity> this[City city]
        {
            get
            {
                return cities[city];
            }
        }
    }
}
