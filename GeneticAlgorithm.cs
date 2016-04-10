using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{

    class GeneticAlgorithm
    {
        public int PopulationSize { get; set; }
        public int TournamentSize { get; set; }
        public double MutationRate { get; set; }
        public int GreedyVariance { get; set; }
        public double GreedyPreferenceRate { get; set; }
        public int MaxGenerations { get; set; }

        public GeneticAlgorithm()
        {
            PopulationSize = 2000;
            TournamentSize = 10;
            MutationRate = 0.015;
            GreedyVariance = 5;
            GreedyPreferenceRate = 0.90;
            MaxGenerations = 100;
        }

        private ProblemAndSolver parent;

        private Random random = new Random();

        public string[] Solve(ProblemAndSolver parent)
        {
            this.parent = parent;
            string[] results = new string[3];
            Dictionary<string, double> improvementRate = new Dictionary<string, double>();

            Population pop = InitPopulation();
            Sequence imdabest = pop.Best;
            Stopwatch timer = new Stopwatch();
            improvementRate[timer.Elapsed.ToString()] = imdabest.Score;
            timer.Start();
            int iterationsSanChangement = 0;
            int count = 0;
            int generations = 0;
            List<Sequence> bestones = new List<Sequence>();
            bestones.Add(imdabest);

            while ((iterationsSanChangement < MaxGenerations && timer.Elapsed.TotalSeconds < 90) ||
                pop.Size < 2)
            {

                generations++;
                CrossPopulation(ref pop);
                MutatePopulation(pop);
                PrunePopulation(pop);
                if (pop.Best.Score < imdabest.Score)
                {
                    iterationsSanChangement = 0;
                    count++;
                    imdabest = pop.Best;
                    improvementRate[timer.Elapsed.ToString()] = imdabest.Score;
                    bestones.Add(imdabest);
                }
                else
                {
                    iterationsSanChangement++;
                }
            }
            timer.Stop();

            parent.BSSF = imdabest.ToRoute();

            Console.WriteLine();

            results[ProblemAndSolver.COST] = imdabest.Score.ToString();    // load results into array here, replacing these dummy values
            results[ProblemAndSolver.TIME] = timer.Elapsed.ToString();
            results[ProblemAndSolver.COUNT] = count.ToString(); //*** should this be 1 or the number of greedy solutions that we find???
            StringBuilder bob = new StringBuilder();
            bob.Append("Time,Best Cost\r\n");
            foreach (var entry in improvementRate)
            {
                bob.Append(entry.Key);
                bob.Append(",");
                bob.Append(entry.Value);
                bob.Append("\r\n");
            }
            bob.Append("Best," + imdabest.Score.ToString() + "\r\n");
            bob.Append("Problem Size," + parent.GetCities().Length + "\r\n");
            bob.Append("Seed," + parent.Seed + "\r\n");
            bob.Append("Generations," + generations.ToString() + "\r\n");
            bob.Append("Total Time," + timer.Elapsed.ToString());
            bob.Append("\r\n");

            System.IO.File.AppendAllText(@"C:\Users\Matt\Desktop\tsp_results.csv", bob.ToString());
            return results;
        }

        private Population InitPopulation()
        {
            Population pop = new Population();

            pop.Add(new Sequence(parent.BSSF));

            City[] cities = parent.GetCities();
            City start = cities[0];
            GreedyMaker maker = new GreedyMaker(cities);
            List<City> greedy = new List<City>();
            City source = start;
            greedy.Add(source);
            while (greedy.Count < cities.Length)
            {
                var destinations = maker[source];
                foreach (var c in destinations)
                {
                    City gcity = c.Destination;
                    if (!greedy.Contains(gcity))
                    {
                        greedy.Add(gcity);
                        source = gcity;
                        break;
                    }
                }
            }
            pop.Add(new Sequence(greedy));


            while (pop.Size < PopulationSize)
            {
                City city = start;
                List<City> gene = new List<City>();
                while (gene.Count < cities.Length)
                {
                    if (!gene.Contains(city)) { gene.Add(city); }
                    city = maker[city][(random.Next() <= GreedyPreferenceRate ?
                        random.Next(GreedyVariance) :
                        random.Next(maker[city].Count))].Destination;
                }
                pop.Add(new Sequence(gene));
            }

            return pop;
        }

        private void CrossPopulation(ref Population pop)
        {
            Population children = new Population();
            children.Add(pop.Best);
            while (children.Size < PopulationSize)
            {

                var mom = TournamentParents(pop);
                var dad = TournamentParents(pop);

                var chillen = Cross(mom, dad);
                children.Add(chillen.Item1);
                children.Add(chillen.Item2);
            }
            pop = children;
        }

        private Sequence TournamentParents(Population pop)
        {
            Population tournament = new Population();
            for (int i = 0; i < TournamentSize; i++)
            {
                tournament.Add(pop.GetCandidate());
            }
            return tournament.Best;
        }

        private Tuple<Sequence, Sequence> Cross(Sequence mom, Sequence dad)
        {
            int swapIndex = random.Next(mom.Length);
            var child1 = new Sequence(mom);
            var child2 = new Sequence(dad);

            var val1 = child1[swapIndex];
            var val2 = child2[swapIndex];

            int index1 = child1.IndexOf(val2);
            int index2 = child2.IndexOf(val1);

            child1[index1] = val1;
            child2[index2] = val2;

            child1.Rescore();
            child2.Rescore();

            var temp = child1[swapIndex];
            child1[swapIndex] = child2[swapIndex];
            child2[swapIndex] = temp;

            return Tuple.Create(child1, child2);
        }

        private void MutatePopulation(Population pop)
        {
            pop.ForEach((seq) => {
                double likeliness = random.NextDouble();
                if (likeliness <= MutationRate)
                {
                    int a = random.Next(seq.Length);
                    int b = random.Next(seq.Length);
                    var t = seq[a];
                    seq[a] = seq[b];
                    seq[b] = t;
                    Console.WriteLine();
                }
                seq.Rescore();
            });
        }

        private void PrunePopulation(Population pop)
        {
            pop.TrimDuplicates();
            //pop.TrimWorst();
            Console.WriteLine();
        }

    }
}