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
        private static ProblemAndSolver parent;

        private static Random random = new Random();

        private class Sequence
        {
            private int[] array;

            private double score;

            public int this[int i]
            {
                get
                {
                    return array[i];
                }
                set
                {
                    array[i] = value;
                }
            }

            public int Length
            {
                get
                {
                    return array.Length;
                }
            }

            public Sequence(int[] thing)
            {
                array = thing;
                score = double.PositiveInfinity;
            }

            public Sequence(Sequence parent)
            {
                array = new int[parent.array.Length];
                for (int i = 0; i < array.Length; i++)
                {
                    array[i] = parent[i];
                }
                score = rescore();
            }

            public Sequence(ProblemAndSolver.TSPSolution best)
            {
                var list = best.Route;
                array = new int[list.Count];
                City[] cities = parent.GetCities();
                int j = 0;
                foreach (City c in list)
                {
                    for (int i = 0; i < cities.Length; i++)
                    {
                        if (cities[i].Equals(c))
                        {
                            array[j] = i;
                        }
                    }
                    j++;
                }
                score = rescore();
            }

            public double Score
            {
                get
                {
                    if (double.IsPositiveInfinity(score))
                    {
                        score = rescore();
                    }
                    return score;
                }
            }

            private double rescore()
            {
                return ToRoute().costOfRoute();
            }

            public int IndexOf(int city)
            {
                for (int i = 0; i < Length; i++)
                {
                    if (array[i] == city) { return i; }
                }
                return -1;
            }

            public ProblemAndSolver.TSPSolution ToRoute()
            {
                ArrayList list = new ArrayList();
                foreach (int i in array)
                {
                    list.Add(parent.GetCities()[i]);
                }
                return new ProblemAndSolver.TSPSolution(list);
            }
        }

        private class Population
        {
            private SortedDictionary<double, List<Sequence>> pop;
            
            public int Size
            {
                get
                {
                    int size = 0;
                    foreach (var d in pop)
                    {
                        size += d.Value.Count;
                    }
                    return size;
                }
            }

            public Population()
            {
                pop = new SortedDictionary<double, List<Sequence>>();
            }

            public void Add(Sequence s)
            {
                if (!pop.ContainsKey(s.Score))
                {
                    pop[s.Score] = new List<Sequence>();
                }
                pop[s.Score].Add(s);
            }

            public Sequence Best
            {
                get
                {
                    if (pop.Count == 0)
                    {
                        return null;
                    }
                    return pop.First().Value.First();
                }
            }

            public List<Sequence> SelectParents()
            {
                List<Sequence> parents = new List<Sequence>();
                Random random = GeneticAlgorithm.random;
                int size = Size;
                List<List<Sequence>> values = Enumerable.ToList(pop.Values);
                for (int i = 0; i < size / 2; i++)
                {
                    var candidates = values[random.Next(values.Count)];
                    var mom = candidates[random.Next(candidates.Count)];
                    parents.Add(mom);
                }
                return parents;
            }

            public void ForEach(Action<Sequence> doit)
            {
                foreach (var entry in pop)
                {
                    foreach (var seq in entry.Value)
                    {
                        doit(seq);
                    }
                }
            }

            public void TrimDuplicates()
            {
                foreach (var entry in pop)
                {
                    if (entry.Value.Count > 1)
                    {
                        var best = entry.Value.First();
                        entry.Value.Clear();
                        entry.Value.Add(best);
                    }
                }
            }

            public void TrimWorst()
            {
                List<List<Sequence>> values = Enumerable.ToList(pop.Values);
                pop.Clear();
                for (int i = 0; i < (3 * values.Count) / 4; i++)
                {
                    pop[values[i].First().Score] = values[i];
                }
            }
        }

        public static string[] Solve(ProblemAndSolver parent)
        {
            GeneticAlgorithm.parent = parent;
            string[] results = new string[3];

            Population pop = InitPopulation();
            Sequence imdabest = pop.Best;
            Stopwatch timer = new Stopwatch();

            timer.Start();
            int iterationsSanChangement = 0;
            int count = 0;
            int generations = 0;
            List<Sequence> bestones = new List<Sequence>();
            bestones.Add(imdabest);
            while ((iterationsSanChangement < 100 && timer.Elapsed.TotalSeconds < 60)||
                pop.Size < 2)
            {

                generations++;
                while (pop.Size < 50)
                {
                    CrossPopulation(pop);
                }
                MutatePopulation(pop);
                PrunePopulation(pop);
                if (pop.Best.Score < imdabest.Score)
                {
                    iterationsSanChangement = 0;
                    count++;
                    imdabest = pop.Best;
                    bestones.Add(imdabest);
                }
                else
                {
                    iterationsSanChangement++;
                }
            }
            timer.Stop();

            parent.BSSF = imdabest.ToRoute();

            Console.WriteLine("Generations:" + generations.ToString());

            results[ProblemAndSolver.COST] = imdabest.Score.ToString();    // load results into array here, replacing these dummy values
            results[ProblemAndSolver.TIME] = timer.Elapsed.ToString();
            results[ProblemAndSolver.COUNT] = count.ToString(); //*** should this be 1 or the number of greedy solutions that we find???

            return results;
        }

        private static Population InitPopulation()
        {
            int pop_size = 50;
            Population pop = new Population();
            int i = 0, swap = -1, temp = -1;
            Random rnd = random;
            City[] Cities = parent.GetCities();

            for (int j = 0; j < pop_size; j++)
            {
                int[] perm = new int[Cities.Length];
                Sequence s = new Sequence(perm);
                do
                {
                    for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                        perm[i] = i;
                    for (i = 0; i < perm.Length; i++)
                    {
                        swap = i;
                        while (swap == i)
                            swap = rnd.Next(0, Cities.Length);
                        temp = perm[i];
                        perm[i] = perm[swap];
                        perm[swap] = temp;
                    }
                    s = new Sequence(perm);
                } while (s.Score == double.PositiveInfinity);
                pop.Add(s);
            }

            return pop;
        }

        private static void CrossPopulation(Population pop)
        {
            City[] cities = parent.GetCities();
            var parents = pop.SelectParents();

            while (parents.Count >= 2)
            {

                var mom = parents[random.Next(parents.Count)];
                var dad = parents[random.Next(parents.Count)];

                var chillen = Cross(mom, dad);

                if (mom.Score < dad.Score)
                {
                    parents.Remove(mom);
                }
                else
                {
                    parents.Remove(dad);
                }

                pop.Add(chillen.Item1);
                pop.Add(chillen.Item2);
                
            }
        }

        private static Tuple<Sequence, Sequence> Cross(Sequence mom, Sequence dad)
        {
            int swapIndex = random.Next(mom.Length);
            var child1 = new Sequence(mom);
            var child2 = new Sequence(dad);

            int val1 = child1[swapIndex];
            int val2 = child2[swapIndex];

            int index1 = child1.IndexOf(val2);
            int index2 = child2.IndexOf(val1);

            child1[index1] = val1;
            child2[index2] = val2;

            int temp = child1[swapIndex];
            child1[swapIndex] = child2[swapIndex];
            child2[swapIndex] = temp;

            return Tuple.Create(child1, child2);
        }
        
        private static double mutation_likelihood = 0.30;

        private static void MutatePopulation(Population pop)
        {
            pop.ForEach((seq) => {
                double likeliness = random.NextDouble();
                if (likeliness <= mutation_likelihood)
                {
                    int a = random.Next(seq.Length);
                    int b = random.Next(seq.Length);
                    int t = seq[a];
                    seq[a] = seq[b];
                    seq[b] = t;
                    Console.WriteLine();
                }
            });
        }

        private static void PrunePopulation(Population pop)
        {
            pop.TrimDuplicates();
            pop.TrimWorst();
            Console.WriteLine();
        }

    }
}
